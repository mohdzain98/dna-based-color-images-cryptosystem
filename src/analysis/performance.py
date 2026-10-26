"""Encryption speed metrics.

Converts measured encryption times into throughput (bytes/sec) and
cycles per byte for a given CPU clock. Used for the time complexity table.
"""

CHANNELS = 3  # RGB


def calculate_et_and_cpb(
    image_dims: list[tuple[int, int]],
    times: list[float],
    cpu_speed_ghz: float = 3.2,
) -> list[dict]:
    """Return encryption throughput and cycles/byte per (height, width) image.

    `times` are encryption times in seconds; 3.2 GHz is the Apple M1 default.
    """
    cpu_speed_hz = cpu_speed_ghz * 1e9
    results = []
    for (height, width), seconds in zip(image_dims, times):
        et_bps = height * width * CHANNELS / seconds
        results.append(
            {
                "Image Size": f"{height}x{width}",
                "Time (s)": seconds,
                "ET (bytes/sec)": round(et_bps, 2),
                "CPB (cycles/byte)": round(cpu_speed_hz / et_bps, 2),
            }
        )
    return results
