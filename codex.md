# Codex Work Summary

## Goal Shift

초기 방향은 "데이터 전처리 함수의 연속"을 recipe로 묶어 재사용 가능한 workflow API로 만드는 것이었다. 이후 범위를 더 명확히 하면서, DataCooker는 정적인 그래프 기반 workflow core library로 두고, StructCooker는 그 위에서 domain-specific pipeline을 구성하는 구조로 정리하기로 했다.

핵심 개념은 다음으로 수렴했다.

- `Read`: 외부 데이터 소스에서 raw input을 가져오는 단계
- `Transform`: recipe/instruction 기반의 정적 변환 단계
- `Write`: 결과를 외부 표현으로 내보내는 단계
- `Runner`: 위 단계를 실제 batch/job/LMDB workflow로 실행하는 계층

이 분리를 둔 이유는, I/O와 순수 변환을 섞으면 재사용성, 테스트 용이성, 배치 실행, 캐시/LMDB 재구성, 실패 복구가 빠르게 복잡해지기 때문이다.

## DataCooker Changes

DataCooker는 core workflow library로 정리하면서 다음 방향으로 다듬었다.

- recipe / loading / executor / runners / CLI / LMDB 기능을 역할별로 분리
- 정적 workflow 실행 API를 library 중심으로 재구성
- LMDB 기능을 `datacooker.lmdb` 아래로 정리
- 공용 DB 관련 유틸은 `datacooker.utils.db` 호환 레이어로 유지
- 문서에서 Read / Transform / Write / Runner 분리를 명시

추가로 StructCooker에서 쓰던 실제 사용 패턴을 맞추기 위해 recipe 입력 타입 처리도 완화했다.

- `str | None` 같은 optional raw annotation 허용
- generic / union annotation 허용 범위 확장
- flattened raw input (`_atom_site.id`)가 있을 때 상위 raw input (`_atom_site`) 요구를 충족하는 방식으로 보완

이 수정으로 StructCooker의 CIF ingest recipe가 DataCooker 위에서 자연스럽게 동작하게 만들었다.

## StructCooker Changes

StructCooker는 DataCooker 위에 올라가는 domain library라는 관점으로 재정리했다.

- 패키지 구조를 `src/structcooker/...` 중심으로 재배치
- `models`는 실제로 biomol 의존 분자 객체이므로 `mols` 관점으로 해석
- constants 성격 파일은 utils 쪽으로 정리
- `mol_type_map.py` 내용은 `utils/mapping.py`로 흡수하는 방향 채택
- instruction 계층을 `instructions/readers`, `instructions/transforms`, `instructions/writers` 식으로 정리
- workflow 파일들은 ingest / analysis 등 의미 기준으로 재배치
- 누락된 데이터 흐름 문서 `docs/data_flow.md` 복구

또한 StructCooker 내부에 남아 있던 공용 DB 제작 코드, 공용 유틸 성격의 코드 중 DataCooker가 가져가야 하는 것들을 지속적으로 분리했다.

## CIF / CCD Verification

StructCooker ingest 정리의 첫 검증 대상으로 CCD LMDB를 선택했다. 이유는 상대적으로 작고 빠르게 검증 가능했기 때문이다.

중간에 `parse_cif.py`가 깨지는 문제가 있었고, 원인은 StructCooker recipe가 쓰는 raw variable annotation (`str | None`)을 DataCooker가 너무 엄격하게 해석하던 점이었다. DataCooker 수정 후에는 `parse_file` 기반 흐름을 유지한 채 정상 동작하도록 복구했다.

그 다음 현재 StructCooker 코드로 CCD 전체를 새 LMDB에 다시 빌드하고, 기존 DB와 비교했다.

비교 결과는 다음과 같았다.

- 기존 DB entries: `49,549`
- 새 DB entries: `49,550`
- 차이 key: 새 DB에만 `UNL` 1건 존재
- 공통 key `49,549`개 전부 deserialize 후 비교: 모두 동일

즉, 현재 CCD ingest는 기존 결과를 사실상 그대로 재현했고, 추가 차이는 `UNL` 1건뿐이다.

## Current Interpretation

현재 상태에서의 해석은 다음과 같다.

- DataCooker는 "정적인 데이터 workflow를 표현하고 실행하는 core library"로 개념이 맞다.
- StructCooker는 그 위에서 biomolecular structure domain recipe / ingest / analysis를 제공하는 library로 두는 것이 맞다.
- I/O를 transform과 완전히 같은 것으로 취급하지 않고 분리한 판단은 유효했다.
- LMDB ingest 기준으로 현재 StructCooker의 주요 파이프라인은 이전 버전과 실질적으로 호환된다.

## Remaining Direction

앞으로의 큰 방향은 다음과 같다.

- DataCooker는 core abstraction의 일관성, 문서, validation, runner ergonomics를 더 다듬는다.
- StructCooker는 recipe 수가 많아질수록 domain별 폴더링, ingest/analysis 목적 분리, instruction naming 정리가 중요하다.
- CCD처럼 작은 ingest부터 old/new 비교를 해가며 다른 workflow도 순차 검증한다.
- `UNL` 1건 차이는 필요하면 별도로 추적해 source data 차이인지 parser 차이인지 확인한다.
