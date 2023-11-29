with distinct_cell_counts as (
  select
    t_inner.query_id,
    ci_inner.cell_type,
    ii_inner.dataset_id,
    min(t_inner.match_score) as min_distance,
    max(t_inner.match_score) as max_distance,
    approx_quantiles(t_inner.match_score, 100)[safe_ordinal(50)] as median_distance,
    avg(t_inner.match_score) as mean_distance,
    count(distinct t_inner.match_cas_cell_index) as count_per_dataset
  from `${temp_table_fqn}` t_inner
  join `${project}.${dataset}.cas_cell_info` ci_inner on ci_inner.cas_cell_index = t_inner.match_cas_cell_index
  join `${project}.${dataset}.cas_ingest_info` ii_inner on ii_inner.cas_ingest_id = ci_inner.cas_ingest_id
  group by t_inner.query_id, ci_inner.cell_type, ii_inner.dataset_id
)
select
  t.query_id,
  ci.cell_type,
  min(t.match_score) as min_distance,
  max(t.match_score) as max_distance,
  approx_quantiles(t.match_score, 100)[safe_ordinal(25)] as p25_distance,
  approx_quantiles(t.match_score, 100)[safe_ordinal(50)] as median_distance,
  approx_quantiles(t.match_score, 100)[safe_ordinal(75)] as p75_distance,
  avg(t.match_score) as mean_distance,
  count(distinct t.match_cas_cell_index) as cell_count,
  (
    select array_agg(
      struct<
      dataset_id string,
      count_per_dataset int,
      min_distance float64,
      max_distance float64,
      median_distance float64,
      mean_distance float64
      > (dataset_id, count_per_dataset, min_distance, max_distance, median_distance, mean_distance)
    )
    from distinct_cell_counts
    where query_id = t.query_id and cell_type = ci.cell_type
  ) as dataset_ids_with_counts
from `${temp_table_fqn}` t
join `${project}.${dataset}.cas_cell_info` ci on t.match_cas_cell_index = ci.cas_cell_index
group by t.query_id, ci.cell_type
order by t.query_id, cell_count desc