select
  t.query_id,
  ci.cell_type,
  min(t.match_score) as min_distance,
  max(t.match_score) as max_distance,
  approx_quantiles(t.match_score, 100)[safe_ordinal(25)] as p25_distance,
  approx_quantiles(t.match_score, 100)[safe_ordinal(50)] as median_distance,
  approx_quantiles(t.match_score, 100)[safe_ordinal(75)] as p75_distance,
  count(*) as cell_count
from `${temp_table_fqn}` t
join `${project}.${dataset}.cas_cell_info` ci on t.match_cas_cell_index = ci.cas_cell_index
group by t.query_id, ci.cell_type
order by t.query_id, cell_count desc