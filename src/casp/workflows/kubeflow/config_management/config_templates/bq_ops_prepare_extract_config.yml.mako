project_id: ${project_id}
obs_columns_to_include: ${obs_columns_to_include}
filters_json_path: ${filters_json_path}
dataset: ${dataset}
extract_table_prefix: ${extract_table_prefix}
fq_allowed_original_feature_ids: ${fq_allowed_original_feature_ids}
extract_bin_size: ${context.get("extract_bin_size", None)}
assign_bin_by_category: ${context.get("assign_bin_by_category", False)}
extract_bin_category_column_name: ${context.get("extract_bin_category_column_name", None)}
bucket_name: ${bucket_name}
extract_bucket_path: ${extract_bucket_path}