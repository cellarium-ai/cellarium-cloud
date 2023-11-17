select test_column,
       test_column_2,
       test_column_3
from `test-project.test_dataset.test_table_name` ttn
join `test-project.test_dataset.test_join_table` tjt on (ttn.id = tjt.id)
where test_column_3 = 5
  and test_column in ('foo',
                      'buz',
                      'ex')
group by 1
limit 10