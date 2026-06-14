# Why Separate Read / Transform / Write?

`DataCooker`에서 `Read / Transform / Write`를 나눈 이유는,
모든 것이 transform이 아니어서가 아니다.

이론적으로는 전부 transform으로 볼 수 있다.

- `file -> object`
- `bytes -> dict`
- `dict -> bytes`
- `object -> object`

하지만 라이브러리 설계에서는 아래 둘을 분리하는 편이 훨씬 낫다.

- `Transform`
  - workflow의 본질적인 계산
  - artifact dependency graph의 내부 의미
- `Read / Write`
  - 외부 세계와 맞닿는 경계
  - 파일 포맷, 압축 방식, DB payload, output materialization

## 왜 필요한가

분리를 안 하면 파이프라인이 보통 다음처럼 불어난다.

- `cif.gz -> parse -> graph -> write_lmdb`
- `cif.gz -> parse -> graph -> write_tsv`
- `lmdb bytes -> deserialize -> graph -> write_lmdb`
- `lmdb bytes -> deserialize -> graph -> write_tsv`

여기서 핵심 계산은 사실 하나다.

- `structure object -> graph feature object`

하지만 입력 표현과 출력 표현이 섞여 있으면,
같은 계산이 `input format x output format` 조합만큼 증식한다.

반대로 분리하면 구조는 이렇게 된다.

- `Read`
  - `cif.gz -> structure object`
  - `lmdb bytes -> structure object`
- `Transform`
  - `structure object -> graph feature object`
- `Write`
  - `graph feature object -> lmdb bytes`
  - `graph feature object -> tsv/json`

그러면:

- 새 입력 포맷이 생기면 `Read`만 추가하면 된다.
- 새 출력 포맷이 생기면 `Write`만 추가하면 된다.
- 핵심 계산은 그대로 재사용된다.

## 한 줄 요약

`Read / Transform / Write` 분리는 예쁘게 나누기 위한 것이 아니라,
**본질적인 계산과 표현 형식/저장 매체/부수효과를 분리해서 변화 전파 범위를 줄이기 위한 설계**다.
