"""Batch processing helpers for static recipe execution over item collections."""

from __future__ import annotations

import logging
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from datacooker._parallel import iter_parallel_chunks
from datacooker.api import parse_dict
from datacooker.protocols import TransformFunc
from datacooker.readers import resolve_key_transform
from datacooker.utils.sharding import resolve_node_config, shard_items

logger = logging.getLogger(__name__)


@dataclass(frozen=True, slots=True)
class BatchProcessResult:
    """Capture the outcome of processing one input mapping."""

    index: int
    input_data: dict[str, Any]
    output: dict[str, Any]
    error: Exception | None = None

    @property
    def succeeded(self) -> bool:
        """Return whether this item completed without an exception."""
        return self.error is None


@dataclass(frozen=True, slots=True)
class BatchProcessReport:
    """Summarize a batch workflow execution."""

    items: tuple[BatchProcessResult, ...]
    assigned: int
    total: int
    node_rank: int
    node_count: int

    @property
    def attempted(self) -> int:
        """Return the number of items processed on this node."""
        return len(self.items)

    @property
    def succeeded(self) -> int:
        """Return the number of successful items."""
        return sum(1 for item in self.items if item.succeeded)

    @property
    def failed(self) -> int:
        """Return the number of failed items."""
        return self.attempted - self.succeeded

    def outputs(self) -> list[dict[str, Any]]:
        """Return the successful output payloads."""
        return [item.output for item in self.items if item.succeeded]

    def errors(self) -> list[Exception]:
        """Return the collected exceptions."""
        return [item.error for item in self.items if item.error is not None]

    def raise_first_error(self) -> None:
        """Raise the first captured exception, if any."""
        errors = self.errors()
        if errors:
            raise errors[0]


def parallel_process(
    data_list: Sequence[Mapping[str, Any]],
    *,
    inputs: Mapping[str, Any] | None,
    recipe: Path,
    key_transform: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    test_run: bool = True,
    node_rank: int | None = None,
    node_count: int | None = None,
    **extra_kwargs: Any,
) -> list[tuple[dict[str, Any], Exception | None]]:
    """Compatibility wrapper returning the historical list-of-tuples shape."""
    report = parallel_process_report(
        data_list,
        inputs=inputs,
        recipe=recipe,
        key_transform=key_transform,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        test_run=test_run,
        node_rank=node_rank,
        node_count=node_count,
        **extra_kwargs,
    )
    return [(item.output, item.error) for item in report.items]


def parallel_process_report(
    data_list: Sequence[Mapping[str, Any]],
    *,
    inputs: Mapping[str, Any] | None,
    recipe: Path,
    key_transform: TransformFunc | None = None,
    chunk_size: int = 10_000,
    n_jobs: int = -1,
    test_run: bool = True,
    node_rank: int | None = None,
    node_count: int | None = None,
    **extra_kwargs: Any,
) -> BatchProcessReport:
    """Parse a collection of input mappings through the same recipe in parallel."""
    if chunk_size < 1:
        msg = f"chunk_size must be >= 1, got {chunk_size}."
        raise ValueError(msg)

    resolved_node_rank, resolved_node_count = resolve_node_config(
        node_rank=node_rank,
        node_count=node_count,
    )
    indexed_data_list = list(enumerate(data_list))
    sharded_items = shard_items(
        indexed_data_list,
        node_rank=resolved_node_rank,
        node_count=resolved_node_count,
    )

    logger.info(
        "Node %d/%d assigned %d/%d entries.",
        resolved_node_rank,
        resolved_node_count,
        len(sharded_items),
        len(data_list),
    )

    if not sharded_items:
        logger.warning(
            "No data assigned to node %d/%d. Exiting early.",
            resolved_node_rank,
            resolved_node_count,
        )
        return BatchProcessReport(
            items=(),
            assigned=0,
            total=len(data_list),
            node_rank=resolved_node_rank,
            node_count=resolved_node_count,
        )

    shared_inputs = dict(inputs or {})
    resolved_key_transform = resolve_key_transform(key_transform=key_transform)

    def _process_item(
        payload: tuple[int, Mapping[str, Any]],
    ) -> BatchProcessResult:
        index, data_dict = payload
        input_copy = dict(data_dict)
        try:
            process_input = dict(input_copy)
            process_input.update(shared_inputs)
            results = parse_dict(
                recipe_path=recipe,
                datadict=process_input,
                key_transform=resolved_key_transform,
                **extra_kwargs,
            )
        except Exception as error:  # noqa: BLE001
            return BatchProcessResult(
                index=index,
                input_data=input_copy,
                output={},
                error=error,
            )
        return BatchProcessResult(
            index=index,
            input_data=input_copy,
            output=results,
        )

    logger.info("To be parsed %d entries.", len(sharded_items))

    if test_run:
        for payload in sharded_items[:10]:
            result = _process_item(payload)
            if result.error is not None:
                raise result.error
        logger.info("Test run successful. Proceeding with full processing...")

    results: list[BatchProcessResult] = []
    for chunk_results in iter_parallel_chunks(
        items=sharded_items,
        chunk_size=chunk_size,
        n_jobs=n_jobs,
        process_item=_process_item,
        prefer="processes",
        log_label="Processing files",
    ):
        results.extend(chunk_results)

    error_count = sum(1 for item in results if item.error is not None)
    logger.info(
        "Processing completed with %d errors out of %d entries.",
        error_count,
        len(results),
    )
    if error_count > 0:
        logger.info("Sample errors:")
        max_logged_errors = 10
        for index, item in enumerate(results):
            if item.error is None:
                continue
            logger.info("Error %d: %s", index + 1, str(item.error))
            if index >= max_logged_errors - 1:
                break

    return BatchProcessReport(
        items=tuple(results),
        assigned=len(sharded_items),
        total=len(data_list),
        node_rank=resolved_node_rank,
        node_count=resolved_node_count,
    )
