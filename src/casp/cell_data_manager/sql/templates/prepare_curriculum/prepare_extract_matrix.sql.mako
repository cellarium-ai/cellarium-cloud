<%!
    from casp.cell_data_manager.sql import mako_helpers
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}__extract_raw_count_matrix`
partition by range_bucket(extract_bin, generate_array(0,${partition_bin_count},${partition_size}))
cluster by extract_bin
as
select ${mako_helpers.parse_column_names(select_columns)},
    array_agg(struct<feature_index int64, raw_counts int64>(ef.cas_feature_index, m.raw_counts)) as feature_data
from `${project}.${dataset}.cas_raw_count_matrix` m
join `${project}.${dataset}.cas_feature_info` fi on (m.cas_feature_index = fi.cas_feature_index)
join `${project}.${dataset}.${extract_table_prefix}__extract_feature_info` ef on (fi.original_feature_id = ef.original_feature_id)
join `${project}.${dataset}.${extract_table_prefix}__extract_cell_info` b on (m.cas_cell_index = b.cas_cell_index)
group by 1,2