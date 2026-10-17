# First offline analysis wiring error, no DFT

The first invocation of `analyze_training_E264.py` stopped before replay with:

```text
ValueError: finite orthonormal free space and unit open branch direction required
```

The caller passed the40-component full-chart soft direction returned by the
training-input loader into the36-component internal optimizer. It now uses
the exact frozen `config['branch_direction_36']`, as the original production
worker does. No dimensional guard, arithmetic tolerance, optimizer, native
data, physical input or calculation was changed. The subsequent ordered
replay consumes all81 records. The actual-material regression exercises this
dimension boundary and would catch the original caller error.
