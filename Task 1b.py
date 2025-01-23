import os
import mysql.connector # Import mysql.connector library to connect and
# interact with MariaDB database
import sqlparse  # Import sqlparse library to parse and split SQL scripts
# into individual statements

# Configuration for MariaDB
db_config = {
    "host": "localhost", # Hostname of database server
    "user": "root", # Username to connect to database
    "password": "password", # Password for database user
    "database": "gene_ontology" # Name of database to use
}

script_directory = os.path.dirname(os.path.abspath(__file__))
task_1_dir = os.path.join(script_directory, "Task 1")
sql_file_path = os.path.join(task_1_dir, "go-dump.sql") # Path and name of
# SQL file


def create_database_if_absent(db_name):
    """
    Create MariaDB database to store SQL file if absent.

    Args:
        db_name (str): The name of the database to create or verify.
    """
    with mysql.connector.connect(
            host=db_config["host"],
            user=db_config["user"],
            password=db_config["password"] # Connect to database server
            # using specified username and password
    ) as conn, conn.cursor() as cursor: # Open database connection and cursor
        cursor.execute(f"CREATE DATABASE IF NOT EXISTS {db_name}") # Create
        # SQL database if not present
        print(f"Database '{db_name}' successfully checked or created.")
        # Print notification of success


def import_sql_file(sql_file):
    """
    Read and execute "go-dump.sql" SQL file.

    Args:
        sql_file (str): Path to SQL file to be imported.
    """
    with mysql.connector.connect(
            host=db_config["host"],
            user=db_config["user"],
            password=db_config["password"],
            database=db_config["database"] # Specify database to use
    ) as conn, conn.cursor() as cursor:
        if not os.path.exists(sql_file): # Check if SQL file exists
            print(f"SQL file not found: {sql_file}") # Notify if file is missing
            return # Exit function

        with open(sql_file, 'r') as file: # Open SQL file for reading ("r")
            sql_script = file.read() # Read file content

        # Parse and split SQL commands using sqlparse
        parsed_statements = sqlparse.split(sql_script) # Split file content
        # into individual SQL statements
        for statement in parsed_statements: # Iterate through SQL statements
            if statement.strip():  # Skip empty statements
                try:
                    cursor.execute(statement) # Execute SQL statement
                except mysql.connector.Error as e: # Handle execution errors
                    print(f"Error executing statement:\n{statement}\n{e}")
                    # Print statement raising execution error
        conn.commit() # Commit changes to database
        print(f"SQL file '{sql_file}' imported successfully.") # Print
        # notification of successful import of SQL file onto MariaDB


def main():
    try:
        # Step 1: Create MariaDB database if absent
        create_database_if_absent(db_config["database"]) # Call
        # create_database_if_absent function

        # Step 2: Import "go-dump.sql" SQL file into database
        import_sql_file(sql_file_path) # Execute SQL script (Call
        # import_sql_file function)
    except mysql.connector.Error as e:
        print(f"Error: {e}") # Print database-related errors


if __name__ == "__main__":
    main()
