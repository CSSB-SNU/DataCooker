# DataCooker V2 Refactor Summary

## Goal

이번 정리의 목표는 `DataCooker`를 정적인 workflow API library로 더 명확하게 재구성하는 것이었다.
핵심 방향은 다음과 같다.

- `Transform`를 recipe graph의 중심 개념으로 유지한다.
- 입출력 경계는 recipe 바깥으로 분리한다.
- LMDB 기능은 generic utility로 유지하되 core workflow API와 같은 언어로 정리한다.
- downstream 프로젝트인 `StructCooker`는 domain layer로 남기고, generic logic은 `DataCooker`로 끌어올린다.

## New V2 Structure

V2 public surface는 네 개의 축으로 정리했다.

- `Read`
  - 입력 파일 로드, payload decode, key path 변환, input normalization
  - module: `datacooker.readers`
- `Transform`
  - static recipe graph execution
  - 기존 `RecipeBook`, `execute`, `parse_file`, `parse_dict` 중심 구조 유지
- `Write`
  - payload encode, result materialization
  - module: `datacooker.writers`
- `Runner`
  - read/transform/write를 합성하는 상위 orchestration
  - module: `datacooker.runners`

## Main API Changes

### Reader side

`ReaderHooks`를 도입했다.

- `loader`
- `adapter`
- `deserializer`
- `key_transform`

추가 helper:

- `dot_path`
- `load_inputs`
- `decode_payload`
- `normalize_input_sequence`
- `resolve_key_transform`

### Writer side

`WriterHooks`를 도입했다.

- `serializer`
- `materializer`

추가 helper:

- `encode_output`
- `write_output`

### Runner side

새 runner API:

- `run_recipe`
- `run_recipe_batch`
- `run_lmdb_extract`

기존 orchestration API:

- `run_workflow`
- `run_parallel_workflow`
- `extract_lmdb_workflow`

위 함수들은 compatibility wrapper로 유지했다.

## Config Migration

config는 flat callable key에서 nested hook 구조를 더 잘 표현하도록 확장했다.

이전:

```yaml
load_func: pkg.module.load_sample
transform_func: pkg.module.dot_transform
project_func: pkg.module.write_result
```

이후:

```yaml
reader:
  loader: pkg.module.load_sample
  key_transform: pkg.module.dot_transform
writer:
  materializer: pkg.module.write_result
```

호환성 때문에 이전 key도 여전히 읽을 수 있게 남겨두었다.

## LMDB Refactor

LMDB pipeline도 같은 vocabulary로 재정리했다.

- `build_lmdb`는 `reader`/`writer` hook을 받을 수 있다.
- `rebuild_lmdb`는 deserialize, adapt, transform, serialize 경계를 분리해서 받는다.
- `extract_lmdb_records`는 LMDB-wide extraction helper로 유지하되 `ReaderHooks` 기반으로 동작한다.
- serializer/deserializer 누락 시 더 명확한 error를 내도록 정리했다.

이 변화로 LMDB 관련 기능이 recipe engine과 별개 규칙으로 움직이지 않고, 같은 conceptual model 위에서 설명된다.

## Public Surface Cleanup

`datacooker.__init__` export를 정리해서 다음 개념을 직접 노출한다.

- protocol alias
  - `FileReader`
  - `DataAdapter`
  - `KeyTransform`
  - `PayloadReader`
  - `PayloadWriter`
  - `ResultWriter`
- hooks
  - `ReaderHooks`
  - `WriterHooks`
- helper API
  - `load_inputs`
  - `decode_payload`
  - `encode_output`
  - `write_output`
- runners
  - `run_recipe`
  - `run_recipe_batch`
  - `run_lmdb_extract`

## StructCooker Coordination

`StructCooker`에서 generic DB/build/extract workflow 쪽은 가능한 한 `DataCooker` vocabulary를 따르도록 정리했다.

확인된 downstream 변경:

- config callable fields를 `reader` / `writer` 구조로 이관
- `scripts/parse_cif.py`의 `load_func` / `transform_func` 사용을 `loader` / `key_transform`으로 변경

의도는 `StructCooker`가 domain recipe와 domain IO에 집중하고, reusable workflow mechanics는 `DataCooker`가 담당하게 만드는 것이다.

## Validation Status

완료:

- `ruff check src tests`

제한:

- `pytest`는 현재 로컬 `pixi` 환경이 깨져 있어 완전 검증하지 못했다.
- 확인된 상태:
  - system python은 `3.9`
  - project requirement는 `>=3.10`
  - `.pixi/envs/default/bin/pytest` shebang이 존재하지 않는 `python3.13` interpreter를 가리킴
  - sandbox 안에서 `pixi run pytest`는 cache/network 제약으로 재생성이 막힘

즉 현재 상태는:

- 정적 lint는 통과
- runtime test는 environment 복구 후 재검증 필요

## Remaining Recommendation

다음 단계는 코드 구조 변경이 아니라 release quality 검증이다.

- `pixi` environment를 복구
- `pytest` 전체 실행
- config examples를 실제 sample workflow로 smoke test
- 필요하면 migration guide를 README에 더 짧게 요약
