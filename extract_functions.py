with open("app.py", encoding="utf-8") as f:
    lines = f.readlines()

def extract_block(start_line_1indexed, out):
    """Extracts a top-level assignment/def block starting at start_line,
    stopping when we hit a line at column 0 that isn't blank, after the
    first line."""
    i = start_line_1indexed - 1
    block = [lines[i]]
    for j in range(i + 1, len(lines)):
        line = lines[j]
        if line.strip() == "":
            block.append(line)
            continue
        if not line[0].isspace():
            break
        block.append(line)
    out.writelines(block)
    out.write("\n" + "=" * 80 + "\n\n")

with open("core_functions_full.txt", "w", encoding="utf-8") as out:
    out.write("=== _KNOWLEDGE_BASE (starts line 1209) ===\n")
    extract_block(1209, out)

    out.write("=== _DEFAULT_PROFILE (starts line 1244) ===\n")
    extract_block(1244, out)

    out.write("=== create_default_microbes (starts line 644) ===\n")
    extract_block(644, out)

print("Done - written to core_functions_full.txt")