import json

if __name__ == "__main__":

    with open("inputs.json", "r") as f:
        inputs_json = json.loads(f.read())

    CONVERT_INPUTS_KEY = "CASIngestData.convert_args"
    old_inputs = inputs_json[CONVERT_INPUTS_KEY]
    new_inputs = []
    c_i = 20000000
    f_i = 20000000

    for _ in range(100):
        for input_obj in old_inputs:
            c_i += 2000000
            f_i += 2000000
            new_inputs.append(
                {"df_filename": input_obj["df_filename"], "cas_cell_index": c_i, "cas_feature_index": f_i}
            )

    inputs_json[CONVERT_INPUTS_KEY] = new_inputs

    with open("repeated_inputs.json", "w+") as f:
        f.write(json.dumps(inputs_json))
