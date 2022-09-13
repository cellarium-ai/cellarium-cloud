PROJECT=gvs-internal
DATASET=cas_phase1_01
PREFIX=demo

# create gene level summary -- scan of full dataset
CREATE OR REPLACE TABLE `gvs-internal.cas_phase1_01.demo__extract_feature_summary`
AS 
SELECT  f.feature_name, 
        sum(m.raw_counts) total_raw_counts,
        COUNT(distinct CASE WHEN m.raw_counts > 0 THEN m.cas_cell_index ELSE null END) cells_with_counts
FROM `gvs-internal.cas_phase1_01.cas_raw_count_matrix` m
JOIN `gvs-internal.cas_phase1_01.cas_feature_info` f ON (m.cas_feature_index = f.cas_feature_index)
GROUP BY f.feature_name


# create subset of features for extract with new local identifiers
# TODO -- join to static list of ensembl ids
CREATE OR REPLACE TABLE `gvs-internal.cas_phase1_01.demo__extract_feature_info`
AS 
SELECT  DENSE_RANK() OVER (ORDER BY s.feature_name ASC) AS cas_feature_index,
        s.feature_name as feature_name,
FROM	`gvs-internal.cas_phase1_01.demo__feature_summary` s 
WHERE s.cells_with_counts >= 3
ORDER BY s.feature_name

#
# TODO: create a cells table with the bin mapping, then just use that mapping in the next query...
CREATE OR REPLACE TABLE `gvs-internal.cas_phase1_01.demo__extract_cell_info`
AS 
SELECT  cas_cell_index,
        CAST(FLOOR(rand() * (select count(1) from `gvs-internal.cas_phase1_01.cas_cell_info`) / 10000) as INT) as extract_bin
FROM	`gvs-internal.cas_phase1_01.cas_cell_info` c
ORDER BY cas_cell_index

# create extract table -- remapping feature identifiers, and including batch identifier
# NOTE: bins are approximately 10k cells, but it's a random sample
CREATE OR REPLACE TABLE `gvs-internal.cas_phase1_01.demo__extract_raw_count_matrix`
PARTITION BY RANGE_BUCKET(id, GENERATE_ARRAY(0,4000,1))
AS 
SELECT  b.extract_bin,
        m.cas_cell_index,
        ef.cas_feature_index,
        m.raw_counts
FROM `gvs-internal.cas_phase1_01.cas_raw_count_matrix` m
JOIN `gvs-internal.cas_phase1_01.cas_feature_info` fi ON (m.cas_feature_index = fi.cas_feature_index)
JOIN `gvs-internal.cas_phase1_01.demo__extract_feature_info` ef ON (fi.feature_name = ef.feature_name)
JOIN `gvs-internal.cas_phase1_01.demo__extract_cell_info` b ON (m.cas_cell_index = b.cas_cell_index)

