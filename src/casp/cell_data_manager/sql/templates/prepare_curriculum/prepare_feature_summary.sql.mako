<%!
    from casp.cell_data_manager.sql import mako_helpers
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}__extract_feature_summary`
as
select f.original_feature_id,
       sum(m.raw_counts) total_raw_counts,
       count(distinct case when m.raw_counts > 0 then m.cas_cell_index else null end) cells_with_counts
from `${project}.${dataset}.cas_raw_count_matrix` m
join `${project}.${dataset}.cas_feature_info` f on (m.cas_feature_index = f.cas_feature_index)
join `${project}.${dataset}.cas_cell_info` c on (m.cas_cell_index = c.cas_cell_index)
% if "dataset_id__eq" in filter_statements.keys() or "dataset_id__in" in filter_statements.keys():
join `${project}.${dataset}.cas_ingest_info` i on (c.cas_ingest_id = i.cas_ingest_id)
% endif
where ${mako_helpers.parse_where_body(filter_statements)}
group by f.original_feature_id