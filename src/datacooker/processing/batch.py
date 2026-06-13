"""Batch processing helpers for static recipe execution over item collections."""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from pathlib import Path
from typing import Any

from datacooker._parallel import iter_parallel_chunks
from datacooker.api import parse_dict
from datacooker.protocols import TransformFunc
from datacooker.utils.sharding import resolve_node_config, shard_items

logger = logging.getLogger(__name__)


def parallel_process(
    data_list: Sequence[Mapping[str, Any]],
    *,
    inputs: Mapping[str, Any] | None,
    recipe: Path,
    transform_func: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    test_run: bool = True,
    node_rank: int | None = None,
    node_count: int | None = None,
    **extra_kwargs: Any,
) -> list[tuple[dict[str, Any], Exception | None]]:
    """Parse a collection of input mappings through the same recipe in parallel."""
    if chunk_size < 1:
        msg = f"chunk_size must be >= 1, got {chunk_size}."
        raise ValueError(msg)

    resolved_node_rank, resolved_node_count = resolve_node_config(
        node_rank=node_rank,
        node_count=node_count,
    )
    sharded_data_list = shard_items(
        list(data_list),
        node_rank=resolved_node_rank,
        node_count=resolved_node_count,
    )

    logger.info(
        "Node %d/%d assigned %d/%d entries.",
        resolved_node_rank,
        resolved_node_count,
        len(sharded_data_list),
        len(data_list),
    )

    if not sharded_data_list:
        logger.warning(
            "No data assigned to node %d/%d. Exiting early.",
            resolved_node_rank,
            resolved_node_count,
        )
        return []

    shared_inputs = dict(inputs or {})

    def _process_item(
        data_dict: Mapping[str, Any],
    ) -> tuple[dict[str, Any], Exception | None]:
        try:
            process_input = dict(data_dict)
            process_input.update(shared_inputs)
            results = parse_dict(
                recipe_path=recipe,
                datadict=process_input,
                transform_func=transform_func,
                **extra_kwargs,
            )
        except Exception as error:  # noqa: BLE001
            return {}, error
        return results, None

    logger.info("To be parsed %d entries.", len(sharded_data_list))

    if test_run:
        for data in sharded_data_list[:10]:
            result = _process_item(data)
            if result[1] is not None:
                raise result[1]
        logger.info("Test run successful. Proceeding with full processing...")

    results: list[tuple[dict[str, Any], Exception | None]] = []
    for chunk_results in iter_parallel_chunks(
        items=sharded_data_list,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_process_item,
        prefer="processes",
        log_label="Processing files",
    ):
        results.extend(chunk_results)

    error_count = sum(1 for _, error in results if error is not None)
    logger.info(
        "Processing completed with %d errors out of %d entries.",
        error_count,
        len(results),
    )
    if error_count > 0:
        logger.info("Sample errors:")
        max_logged_errors = 10
        for index, (_, error) in enumerate(results):
            if error is None:
                continue
            logger.info("Error %d: %s", index + 1, str(error))
            if index >= max_logged_errors - 1:
                break

    return results
