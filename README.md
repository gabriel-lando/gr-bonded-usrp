# gr-bonded-usrp

GNU Radio OOT module that exposes the UHD bonded C++ receiver engine
(`uhd::usrp::bonded::bonded_receiver`) as a source block.

## Status

Initial implementation in progress.

Implemented:

- `bonded_source` C++ block wrapping `bonded_receiver`
- Python bindings (`gnuradio.bonded_usrp`)
- GRC block definition
- Optional per-device frequency plan (`freq_plan_csv`)

Current assumptions:

- Output ports are fixed to 4 (2 devices x 2 channels)
- Existing WiFi flow uses port 0 (dev0ch0) and port 2 (dev1ch0)

## Build

```bash
cd /home/gabriel/uhd/host/gr-bonded-usrp
mkdir -p build && cd build
cmake ..
make -j$(nproc)
sudo make install
sudo ldconfig
```

## Python usage

```python
from gnuradio import bonded_usrp

src = bonded_usrp.bonded_source(
    "30B56D6,30DBC3C",   # serials_csv
    10e6,                 # per-device rate
    100e6,                # uniform freq (ignored when freq_plan_csv is set)
    40.0,                 # gain dB
    "external",          # clock_source
    "external",          # time_source
    False,                # strict_lock
    5.0,                  # lock_timeout
    1.0,                  # fetch_timeout
    "A:B",               # subdev
    "95e6,105e6",        # freq_plan_csv (optional)
)
```
