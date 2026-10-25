# forcing

The forcing component exports one 0-D coefficient per name of `[forcing] coefs`
(contracts-design section 2.7); the names and units come from the par file, so
there is no row to list. A 0-D export without a row is written lossless by the
default writer; tsgen writes its own series (`coef_<name>`, rates) in stream 1D.
See kryos-variables-ice.md for the columns.

| id | variable | dimensions | units | long_name | prec |
|----|----------|------------|-------|-----------|------|
