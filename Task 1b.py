import os
import mysql.connector
import sqlparse  # For proper SQL parsing

# Configuration for MariaDB
db_config = {
    "host": "localhost",
    "user": "root",
    "password": "Diegito16599!",
    "database": "gene_ontology"
}

# Path to SQL file in the script's directory
script_directory = os.path.dirname(os.path.abspath(__file__))
task_1_dir = os.path.join(script_directory, "Task 1")
sql_file_path = os.path.join(task_1_dir, "go-dump.sql")


def create_database_if_not_exists(db_name):
    """Create the database if it doesn't already exist."""
    with mysql.connector.connect(
            host=db_config["host"],
            user=db_config["user"],
            password=db_config["password"]
    ) as conn, conn.cursor() as cursor:
        cursor.execute(f"CREATE DATABASE IF NOT EXISTS {db_name}")
        print(f"Database '{db_name}' checked or created.")


def import_sql_file(sql_file):
    """Read and execute the SQL file."""
    with mysql.connector.connect(
            host=db_config["host"],
            user=db_config["user"],
            password=db_config["password"],
            database=db_config["database"]
    ) as conn, conn.cursor() as cursor:
        if not os.path.exists(sql_file):
            print(f"SQL file not found: {sql_file}")
            return

        with open(sql_file, 'r') as file:
            sql_script = file.read()

        # Parse and split SQL commands using sqlparse
        parsed_statements = sqlparse.split(sql_script)
        for statement in parsed_statements:
            if statement.strip():  # Skip empty statements
                try:
                    cursor.execute(statement)
                except mysql.connector.Error as e:
                    print(f"Error executing statement:\n{statement}\n{e}")
        conn.commit()
        print(f"SQL file '{sql_file}' imported successfully.")


def main():
    try:
        # Create the database if it doesn't exist
        create_database_if_not_exists(db_config["database"])

        # Import the SQL file into the database
        import_sql_file(sql_file_path)
    except mysql.connector.Error as e:
        print(f"Error: {e}")


if __name__ == "__main__":
    main()
