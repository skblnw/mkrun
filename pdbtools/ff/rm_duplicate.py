import argparse
import os

def remove_duplicate_lines_in_sections(input_file, output_file, rough_duplicate_search=False):
    with open(output_file, 'w') as output_f:
        with open(input_file, 'r') as input_f:
            current_section = None
            seen_lines = {}
            removed_lines_count = 0
            num_fields_to_check = None

            for line in input_f:
                if line.startswith("["):
                    if current_section:
                        print(f"Section '{current_section}': Removed {removed_lines_count} duplicate lines.")
                        removed_lines_count = 0

                    current_section = line.strip()
                    seen_lines[current_section] = set()
                    output_f.write(line)
                elif current_section:
                    clean_line = line.strip()
                    if clean_line.startswith(";"):  # Skip comment lines
                        output_f.write(line)
                        continue

                    if num_fields_to_check is None:  # Determine N+1 for the section
                        fields = clean_line.split()
                        for field in fields:
                            if field.isdigit():
                                num_fields_to_check = int(field) + 1
                                break

                    if rough_duplicate_search:
                        fields = clean_line.split(None, num_fields_to_check)
                        fields_key = " ".join(fields[1:num_fields_to_check])
                        if fields_key not in seen_lines[current_section]:
                            seen_lines[current_section].add(fields_key)
                            output_f.write(line)
                        else:
                            removed_lines_count += 1
                    else:
                        if clean_line not in seen_lines[current_section]:
                            seen_lines[current_section].add(clean_line)
                            output_f.write(line)
                        else:
                            removed_lines_count += 1

            if current_section:
                print(f"Section '{current_section}': Removed {removed_lines_count} duplicate lines.")

def main():
    parser = argparse.ArgumentParser(description="Remove duplicate lines within sections.")
    parser.add_argument("input_file", help="Path to the input text file")
    parser.add_argument("-o", "--output_file", help="Path to the output text file", default=None)
    parser.add_argument("-r", "--rough_duplicate_search", help="Perform rough duplicate search based on the first N fields", action="store_true")
    args = parser.parse_args()

    if args.output_file is None:
        output_filename = "new_" + os.path.basename(args.input_file)
        args.output_file = os.path.join(os.path.dirname(args.input_file), output_filename)

    remove_duplicate_lines_in_sections(args.input_file, args.output_file, args.rough_duplicate_search)

if __name__ == "__main__":
    main()
