# Source hashes and line endings

> **Historical record.** These runs were made with the sources of commit `1d7c1f8`. Five of
> the files below were changed on 10 November 2026; the current sources and runs are
> documented in [`../rerun-2026-11-10/SOURCE_HASHES.md`](../rerun-2026-11-10/SOURCE_HASHES.md).
> "Committed file" below means the file at commit `1d7c1f8`.

`reproduction.json` and the certificate JSON files record SHA-256 hashes of the source files
as they were on disk when the runs were made, on Windows. Four files then had CRLF line
endings. Git stores every file with LF endings, so a checkout may give LF or CRLF depending on
`core.autocrlf`. The code is the same either way. The table lists the hash of each committed
file under both conventions; the recorded hash is the one marked.

| file | recorded with | SHA-256 with LF endings | SHA-256 with CRLF endings |
|---|---|---|---|
| `certify.py` | CRLF | `92a35cc5a1d11de9f40fdae224de04706f5ef74ef8ecd328a1c89e07ba6b3ad1` | `c4da8dccb8a3736759c81a68364489a4bbdbfe4c56ee2fac06481786350cfd6c` |
| `base_triangle.py` | LF | `44acee044769a6d0e3b1bb949995e2156a64058e8d6164576858d048f13719aa` | `c14fd9d90e3045769a15fb6ea7949af4748c05c866f4a744d27653b31f3ba2a4` |
| `verify.py` | CRLF | `99d2e07bd9cbf095dde79ecf4330b6e62f49af5f9194b12f971324f3dcd0cd52` | `569581f5fadb972a92fbb49deeb063d19c466fdf98d24d80a6a319a830dbbd21` |
| `interval.py` | CRLF | `eaffdf8fc2de71f933dbc17cf747c3dccd1d8b136a394dab306f3bd0db99e620` | `8c3e4fdf7051477fc621b2a8cdae487b3844b7bfbf4d0f35baf4fbc9d2d9d317` |
| `ad.py` | CRLF | `3dadf28b046ca08b79519b692dab54619e1394da0ce18210341d2d078ebf0e71` | `03964d91f2020200b82f7f5f9304a4ab3021c555557535df1fae680c937f6031` |
| `domains.py` | LF | `7351f66697a31baf2add1e705d9e68dc4cfb694c66bbc50c107514fff3f57338` | `efc39b81f9177fe324405869b069ed06c5904847a0696d804c9c0c75c31e7313` |
| `provenance.py` | LF | `8b20470f4d13a7dd52bc82765b1a5b54efe067908520e71e7595db0b4d428434` | `54782009690042098bda10e17ab2ddce397a64c25e64382ae906702bca83abd9` |
| `closedform.py` | LF | `dbd35c3b221587e36148a905820c2351b062d0527173ac41e275ceecb36b9ac3` | `f606376265f9eb4067493b9bd871bed4d918c1c5a6fbfd51b0c4627090932cc8` |
| `reproduce.py` | LF | `beaf09ec2d4f41baf6435de448c5d3a57b4ad61940aa3b7059d445f637d70c7f` | `21b26b2684173d6608f65d1822131b1ef98d8f86ccadaf7d11835cde3f33b494` |
| `tests/test_enclosures.py` | LF | `d50d6d791168ef6297176547338c30ba19449f46623a952400fefd47e97d2c49` | `cf2f5d4ead4e915ba3969ac48caae179ba9a1abc5e31a84fc83e9224f2bd832c` |
| `tests/test_rigorous_arithmetic.py` | LF | `da85469d314dadf7657fdf0d228cbe4f031bdc1198d5121d5c99f3bc213832fd` | `f97df3e5ae00535ba717e70bd86cc74e878e1a5ff44aca356de6117fa70eaab3` |

To check a checkout independently of its line endings:

```python
import hashlib
import pathlib

def h(path, crlf=False):
    data = pathlib.Path(path).read_bytes().replace(b'\r\n', b'\n')
    return hashlib.sha256(data.replace(b'\n', b'\r\n') if crlf else data).hexdigest()

print(h('certify.py', crlf=True))  # compare with the table above
```
