import json

if __name__ == "__main__":
    with open("inputs_extract.json", "r") as f:
        inputs = json.loads(f.read())
        bin_ranges = []
        start_bin, end_bin = 0, 33118
        step_size = 96
        j, k = 0, 95

        for i in range(end_bin // step_size):
            bin_ranges.append([j, k])
            j += step_size
            k += step_size

        bin_ranges.append([j, end_bin])
        inputs["CASExtractBQ.bin_borders"] = bin_ranges

        with open("inputs_benchmark_400m_2.json", "w") as f:
            f.write(json.dumps(inputs))
