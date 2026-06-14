"""Public core API for DataCooker's static workflow engine."""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError, version

try:
    __version__ = version(__name__)
except PackageNotFoundError:
    __version__ = "unknown"

from .api import (
    describe,
    execute,
    parse_dict,
    parse_file,
    visualize,
)
from .cache import ExecutionContext
from .errors import (
    CycleError,
    DataCookerError,
    DuplicateTargetError,
    InstructionOutputError,
    InvalidRecipeError,
    MissingDependencyError,
    MissingTargetError,
    StepExecutionError,
    UnknownTargetError,
)
from .executor import Cooker
from .loading import load_recipe
from .processing import BatchProcessReport, BatchProcessResult, parallel_process_report
from .protocols import (
    ConvertFunc,
    DeserializeFunc,
    KeyFunc,
    LoadFunc,
    ProjectFunc,
    SerializeFunc,
    TransformFunc,
)
from .readers import (
    InputHooks,
    PayloadHooks,
    ReaderHooks,
    decode_payload,
    dot_path,
    load_inputs,
)
from .recipe import Inputs, Recipe, RecipeBook, RecipeError, Variable, variable
from .runners import run_lmdb_extract, run_recipe, run_recipe_batch
from .transforms import (
    extract_float_single,
    key_stack,
    merge_dict,
    single_value_instruction,
)
from .utils.importing import resolve_object
from .utils.paths import scan_paths
from .utils.sharding import resolve_node_config, shard_items
from .writers import WriterHooks, encode_output, write_output

__all__ = [
    "BatchProcessReport",
    "BatchProcessResult",
    "ConvertFunc",
    "Cooker",
    "CycleError",
    "DataCookerError",
    "DeserializeFunc",
    "DuplicateTargetError",
    "ExecutionContext",
    "InputHooks",
    "Inputs",
    "InstructionOutputError",
    "InvalidRecipeError",
    "KeyFunc",
    "LoadFunc",
    "MissingDependencyError",
    "MissingTargetError",
    "PayloadHooks",
    "ProjectFunc",
    "ReaderHooks",
    "Recipe",
    "RecipeBook",
    "RecipeError",
    "SerializeFunc",
    "StepExecutionError",
    "TransformFunc",
    "UnknownTargetError",
    "Variable",
    "WriterHooks",
    "__version__",
    "decode_payload",
    "describe",
    "dot_path",
    "encode_output",
    "execute",
    "extract_float_single",
    "key_stack",
    "load_inputs",
    "load_recipe",
    "merge_dict",
    "parallel_process_report",
    "parse_dict",
    "parse_file",
    "resolve_node_config",
    "resolve_object",
    "run_lmdb_extract",
    "run_recipe",
    "run_recipe_batch",
    "scan_paths",
    "shard_items",
    "single_value_instruction",
    "variable",
    "visualize",
    "write_output",
]
