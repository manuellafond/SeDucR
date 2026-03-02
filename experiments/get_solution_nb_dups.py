import pandas as pd
import sys


#WARNING: this is all chatgpt

def main(csv_file):
    # Read the CSV file
    df = pd.read_csv(csv_file)

    # Group by 'method' and 'duprate' and compute mean and max of 'solution_nbdups'
    grouped = (
        df.groupby(['method', 'duprate'])['solution_nbdups']
          .agg(['mean', 'max'])
          .reset_index()
    )

    # Rename columns for clarity (optional)
    grouped = grouped.rename(columns={
        'mean': 'avg_solution_nbdups',
        'max': 'max_solution_nbdups'
    })

    # Print results
    print(grouped)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python script.py <input.csv>")
        sys.exit(1)

    main(sys.argv[1])