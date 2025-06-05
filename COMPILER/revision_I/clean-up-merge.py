#!/usr/bin/env python3
import os

def remove_last_two_lines(lines, comparison_strings):
  if len(lines) >= 2 and lines[-2].strip() in comparison_strings and lines[-1].strip().upper() == "END":
    removed = lines[-2].strip()
    lines = lines[:-2]
  
    print(f"Last two lines matching one of {removed} removed")
  else:
    print(f"No matching lines found at the end")
  return lines

file_path = os.path.join('Input_total_commented','C3.dic')
print("Reading '" + file_path + "'")
with open(file_path, 'r', encoding='cp1252') as file:
  lines = file.readlines()

print("Removing unused lines")
comparison_strings = ['ABSTRACTORS', 'CORRECTIONS', 'KEQ_REF_ABSTRACTIONS']
for cmp_str in comparison_strings:
  lines = remove_last_two_lines(lines, comparison_strings)

print("Writing '" + file_path + "'")
with open(file_path, 'w', encoding='cp1252') as file:
  file.writelines(lines)

