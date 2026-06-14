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
    parse,
    parse_dict,
    parse_file,
    rebuild,
    visualize,
)
from .cache import ExecutionContext, ParsingCache
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
from .orchestration import (
    extract_lmdb_workflow,
    run_parallel_workflow,
    run_workflow,
)
from .processing import BatchProcessReport, BatchProcessResult, parallel_process_report
from .protocols import (
    ConvertFunc,
    DataAdapter,
    DeserializeFunc,
    FileReader,
    KeyFunc,
    KeyTransform,
    LoadFunc,
    PayloadReader,
    PayloadWriter,
    ProjectFunc,
    ResultWriter,
    SerializeFunc,
    TransformFunc,
)
from .readers import ReaderHooks, decode_payload, dot_path, load_inputs
from .recipe import Inputs, Recipe, RecipeBook, RecipeError, Variable, variable
from .runners import run_lmdb_extract, run_recipe, run_recipe_batch
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
    "DataAdapter",
    "DataCookerError",
    "DeserializeFunc",
    "DuplicateTargetError",
    "ExecutionContext",
    "FileReader",
    "Inputs",
    "InstructionOutputError",
    "InvalidRecipeError",
    "KeyFunc",
    "KeyTransform",
    "LoadFunc",
    "MissingDependencyError",
    "MissingTargetError",
    "ParsingCache",
    "PayloadReader",
    "PayloadWriter",
    "ProjectFunc",
    "ReaderHooks",
    "Recipe",
    "RecipeBook",
    "RecipeError",
    "ResultWriter",
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
    "extract_lmdb_workflow",
    "load_inputs",
    "load_recipe",
    "parallel_process_report",
    "parse",
    "parse_dict",
    "parse_file",
    "rebuild",
    "resolve_node_config",
    "resolve_object",
    "run_lmdb_extract",
    "run_parallel_workflow",
    "run_recipe",
    "run_recipe_batch",
    "run_workflow",
    "scan_paths",
    "shard_items",
    "variable",
    "visualize",
    "write_output",
]
