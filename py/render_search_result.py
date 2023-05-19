import pandas as pd
import yaml

# Load config file
with open("config/render_result.yaml", "r") as f:
    config = yaml.safe_load(f)


# Define a function to generate HTML links
def generate_link(row):
    link = f'<a href="{row[config["linked_column"]]}.html">{row[config["linked_column"]]}</a>'
    return link


def order_df_columns(df):
    # Get the column names and pop the last column
    columns = df.columns.tolist()
    last_column = columns.pop()
    # Reorder the column names and insert the last column at the beginning
    new_columns = [last_column] + columns
    # Reindex the DataFrame with the new column order
    df_out = df.reindex(columns=new_columns)

    return df_out


def run(blastout_dict):
    # Import blastout parsed dictionary into pandas dataframe
    df = pd.DataFrame.from_dict(blastout_dict[list(blastout_dict.keys())[0]], orient='index')

    # Create a column for hit IDs
    df[config["linked_column"]] = df.index
    # Reorder the dataframe so the hit IDs are placed first
    df = order_df_columns(df)
    # Apply the function to links to the hit ID columns
    df[config["linked_column"]] = df.apply(generate_link, axis=1)

    # Convert the DataFrame to HTML
    html = df.to_html(escape=False, index=False)

    # Save the HTML to a file
    with open('output.html', 'w') as f:
        f.write(html)
