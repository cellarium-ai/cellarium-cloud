<%
    from casp.scripts.bq_ops import constants

    feature_biotype = constants.MRNA_FEATURE_BIOTYPE_NAME
%>
update `${project}.${dataset}.cas_cell_info` ci
set ci.total_mrna_umis = m.raw_counts_total
from (
    select cas_cell_index, sum(raw_counts) as raw_counts_total
    from `${project}.${dataset}.cas_raw_count_matrix` m
    join `${project}.${dataset}.cas_feature_info` fi on (fi.cas_feature_index = m.cas_feature_index)
    where fi.feature_biotype = '${feature_biotype}'
    group by cas_cell_index
) as m
where ci.cas_cell_index = m.cas_cell_index
