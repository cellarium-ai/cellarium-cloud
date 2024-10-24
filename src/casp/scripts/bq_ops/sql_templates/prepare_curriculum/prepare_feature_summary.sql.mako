<%!
    from casp.data_manager.sql import mako_helpers as mh
%>
create or replace table `${project}.${dataset}.${extract_table_prefix}__extract_feature_summary`
as
select f.original_feature_id,
       sum(m.raw_counts) as total_raw_counts,
       count(distinct case when m.raw_counts > 0 then m.cas_cell_index else null end) as cells_with_counts
from `${project}.${dataset}.cas_raw_count_matrix` m
join `${project}.${dataset}.cas_feature_info` f on (m.cas_feature_index = f.cas_feature_index)
join `${project}.${dataset}.cas_cell_info` c on (m.cas_cell_index = c.cas_cell_index)
% if "dataset_id" in " ".join(filter_statements.keys()):
join `${project}.${dataset}.cas_ingest_info` i on (c.cas_ingest_id = i.cas_ingest_id)
% endif
${mh.where(filter_statements)}
group by f.original_feature_id