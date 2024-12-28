#!/usr/bin/env python3

import math
from pathlib import Path
from typing import Tuple, List

def read_total_line(filepath: str) -> str:
    """
    Read the line containing 'TOTAL:' from a file.
    
    Args:
        filepath (str): Path to the input file
        
    Returns:
        str: The line containing 'TOTAL:'
        
    Raises:
        FileNotFoundError: If the input file doesn't exist
        ValueError: If no 'TOTAL:' line is found
    """
    try:
        with open(filepath, 'r') as f:
            for line in f:
                if "TOTAL:" in line:
                    return line
        raise ValueError(f"No 'TOTAL:' line found in {filepath}")
    except FileNotFoundError:
        raise FileNotFoundError(f"File not found: {filepath}")

def extract_values(line: str) -> Tuple[float, float]:
    """
    Extract value and deviation from a result line.
    
    Args:
        line (str): Input line containing the values
        
    Returns:
        Tuple[float, float]: (value, deviation)
        
    Raises:
        ValueError: If the line format is invalid
    """
    try:
        values: List[float] = [float(value) for value in line.split()[16:20:2]]
        return values[0], values[1]
    except (IndexError, ValueError):
        raise ValueError("Invalid line format")

def calculate_difference(val1: float, dev1: float, val2: float, dev2: float) -> Tuple[float, float]:
    """
    Calculate the difference and combined error between two values with their deviations.
    
    Args:
        val1 (float): First value
        dev1 (float): First deviation
        val2 (float): Second value
        dev2 (float): Second deviation
        
    Returns:
        Tuple[float, float]: (difference, combined error)
    """
    difference = val1 - val2
    combined_error = math.sqrt(dev1**2 + dev2**2)
    return difference, combined_error

def main():
    # Define file paths
    paths = {
        'free': Path("free/xvg/results.txt"),
        'complex': Path("bound/xvg/results.txt"),
        'output': Path("combined_results.txt")
    }
    
    try:
        # Read and process input files
        complex_line = read_total_line(str(paths['complex']))
        free_line = read_total_line(str(paths['free']))
        
        # Extract values
        complex_val, complex_dev = extract_values(complex_line)
        free_val, free_dev = extract_values(free_line)
        
        # Calculate result
        result, error = calculate_difference(complex_val, complex_dev, 
                                          free_val, free_dev)
        
        # Format output
        formatted_result = f"{result:7.2f} +- {error:4.2f}\n"
        
        # Write results
        paths['output'].parent.mkdir(parents=True, exist_ok=True)
        with open(paths['output'], "w") as f:
            f.write(formatted_result)
            
    except (FileNotFoundError, ValueError) as e:
        print(f"Error: {e}")
        exit(1)

if __name__ == "__main__":
    main()
    