# BQ Ops Services
These are wrappers around the code in `casp.bq_scripts`. These wrappers try to follow the conventions of the old scripts from `casp.bq_scripts` to avoid massive code base change in one time.   Each of the `bq_ops` service has it's own cromwell workflow in `casp.cromwell.bq_ops`.

## Test Data Cycle
A script that tests data infrastructure if it could produce relevant chunks after harmonization and randomization. 
It is excluded from normal tests as it is not a unit test. Just run this test manually:
```
$ python casp/services/bq_ops/test_data_cycle/main.py
```
You can run this as a widdle workflow. Test raises error after all test cases are finished. This would insure you widdle workflow would become marked as `Failed`