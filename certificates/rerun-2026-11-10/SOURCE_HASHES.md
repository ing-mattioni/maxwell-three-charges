# Source hashes and line endings

`reproduction.json`, the certificate JSON files in this directory and the JSON files in
`independent/` record SHA-256 hashes of the source files as they were on disk when the runs
were made, on Windows, on 10 November 2026. Git stores every file with LF endings, so a
checkout may give LF or CRLF depending on `core.autocrlf`. The code is the same either way.
The table lists the hash of each committed file under both conventions; the recorded hash is
the one marked.

Compared with `../review-2026-11-10/`, three files have new content: `interval.py`,
`tests/test_rigorous_arithmetic.py` and `independent/evaluator.py`, which now check the
arithmetic environment. `independent/certify_domain.py` and `independent/certify_base.py`
have the same content as before; their recorded hashes differ only because the files were
on disk with CRLF endings this time.

| file | recorded with | SHA-256 with LF endings | SHA-256 with CRLF endings |
|---|---|---|---|
| `certify.py` | CRLF | `054bfc7211c321ce00e060cbc3c35d17a7215d812a93326be80b62872333b665` | `e8556f7b031f5bb50fe307fa0ddddf5a35f8baeeeccf4886a932c52efa880686` |
| `base_triangle.py` | LF | `44acee044769a6d0e3b1bb949995e2156a64058e8d6164576858d048f13719aa` | `c14fd9d90e3045769a15fb6ea7949af4748c05c866f4a744d27653b31f3ba2a4` |
| `verify.py` | CRLF | `58d8140ec6b76fd02aa5cb3268e3e98db1a3e983f1e9f23b3a29c208735258ad` | `1f88503ae63c8f0fc2a407cd970cf6a6f58d033ef9ffdfcf0c2cda8d8fa2653d` |
| `interval.py` | CRLF | `bac648f131b741ab0b3db85cf30f98f3d11c1496275b07a0de9cbabdefafb17c` | `6f429cc05393808a2d4f577b55d5fa652971cb455222f68ef0608e89bf72d4b7` |
| `ad.py` | CRLF | `9a2c80d16815f867b8d858cac28277cfbb5eae38983d457dc9de121ec5a2229a` | `a84c45b4f5da5d2de51786afa21e4d709c94913bcb3f50793ca096d0f65c879c` |
| `domains.py` | LF | `7351f66697a31baf2add1e705d9e68dc4cfb694c66bbc50c107514fff3f57338` | `efc39b81f9177fe324405869b069ed06c5904847a0696d804c9c0c75c31e7313` |
| `provenance.py` | LF | `8b20470f4d13a7dd52bc82765b1a5b54efe067908520e71e7595db0b4d428434` | `54782009690042098bda10e17ab2ddce397a64c25e64382ae906702bca83abd9` |
| `closedform.py` | LF | `dbd35c3b221587e36148a905820c2351b062d0527173ac41e275ceecb36b9ac3` | `f606376265f9eb4067493b9bd871bed4d918c1c5a6fbfd51b0c4627090932cc8` |
| `reproduce.py` | LF | `beaf09ec2d4f41baf6435de448c5d3a57b4ad61940aa3b7059d445f637d70c7f` | `21b26b2684173d6608f65d1822131b1ef98d8f86ccadaf7d11835cde3f33b494` |
| `tests/test_enclosures.py` | LF | `d50d6d791168ef6297176547338c30ba19449f46623a952400fefd47e97d2c49` | `cf2f5d4ead4e915ba3969ac48caae179ba9a1abc5e31a84fc83e9224f2bd832c` |
| `tests/test_rigorous_arithmetic.py` | CRLF | `11bd0e70c18682bdae088f51bf5f5cb22994083151cb465118b10358a3b87517` | `8e51a69aca2a2d546addc32d576bc054988deb36d31d60e3ef9f7492df0e621f` |
| `independent/certify_base.py` | CRLF | `1ba5bc1344ce77595e720a4a188172317d72f9952dce1d3f34cff940fdd67f85` | `e48301a646295335211ed4ef6242e23e3df795b0ccfd8222de560712ffd8def6` |
| `independent/certify_domain.py` | CRLF | `70e4e3ef5ffcdb52bd262dbeed35414ce942bcd7829be02126fbbfebdd3ed5fb` | `50d1740f1c9daa2e6ef3a596d5722669b80a04d4e4ec6f78df8f9f28f2b20240` |
| `independent/evaluator.py` | CRLF | `a2c7505257e071b253fd88c7b171a4d0e7fe92968e7d3a75a7c6757d180a84ff` | `a3f001354486da997fdaf22571ff876d70239f923a61a04e636f6d075cc60a19` |

To check a checkout independently of its line endings:

```python
import hashlib
import pathlib

def h(path, crlf=False):
    data = pathlib.Path(path).read_bytes().replace(b'\r\n', b'\n')
    return hashlib.sha256(data.replace(b'\n', b'\r\n') if crlf else data).hexdigest()

print(h('certify.py', crlf=True))  # compare with the table above
```
