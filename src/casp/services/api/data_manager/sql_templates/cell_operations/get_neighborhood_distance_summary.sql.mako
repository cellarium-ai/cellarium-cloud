<%
"""
Description:
  This template generates a SQL query that returns a summary of the distances for each neighbor for the query cells.

Parameters:
  - temp_table_fqn: Fully qualified name of the temporary table used in the query.
  - project: Google Cloud project name.
  - dataset: Name of the dataset within the BigQuery project.

Tables Accessed:
  - Temporary table as defined by `temp_table_fqn`.
  - cas_cell_info in the specified project and dataset.

Output Format:
  The query returns the following columns:
  - query_id: Identifier of the query.
  - cell_type: Type of the cell.
  - min_distance: Minimum match score.
  - max_distance: Maximum match score.
  - p25_distance: 25th percentile of match scores.
  - median_distance: Median of match scores.
  - p75_distance: 75th percentile of match scores.
  - cell_count: Count of cells.

Notes:
  - The match scores are calculated and grouped by query_id and cell_type.
  - The query uses APPROX_QUANTILES for efficiency in calculating percentile values.
"""
%>
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
