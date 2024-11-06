# UCF2SBOL

This README provides a detailed guide to installing, configuring, and using UCF2SBOL.

UCF2SBOL is a Java-based tool designed to convert User Constraint File (UCF) data into the SBOL (Synthetic Biology Open Language) format and then upload the SBOL files to Synbiohub.

This tool facilitates the standardized representation of biological parts and circuits in a format that can be stored and shared through SynBioHub.

## Table of Contents

1. [Features](#features)

2. [Requirements](#requirements)

3. [Installation](#installation)

4. [Usage](#usage)

5. [Command-Line Arguments](#command-line-arguments)

## Features

- Converts UCF JSON files to SBOL format.

- Supports upload to SynBioHub for sharing and storing SBOL documents.

- Customizable collection details (e.g., ID, version, name, description, PubMed ID).

## Requirements

- **Java Development Kit** (https://www.oracle.com/java/technologies/downloads/)

If you are already running Eclipse or Intellij IDEA, you should have already JDK installed.

- **libSBOLj** library: `libSBOLj-2.4.0-withDependencies.jar`

Download from here: https://github.com/SynBioDex/libSBOLj/releases/tag/v2.4.0

## Installation

1. **Clone or Download the UCF2SBOL Code**:
   
   - git clone https://github.com/MyersResearchGroup/UCF2SBOL.git

2. **Open the Project**:
- Open IntelliJ IDEA.

- Select **File** > **Open** and navigate to the directory where you saved the UCF2SBOL code.

#### Adding External Libraries

To run UCF2SBOL, you need to add the `libSBOLj-2.4.0-withDependencies.jar` library.

1. **Open Module Settings**:
- Right-click your project in the **Project** view and select **Open Module Settings** (or press `F4`).
2. **Go to Libraries**:
- In the **Project Structure** window, select **Libraries** from the left sidebar.
3. **Add the JAR File**:
- Click the **+** icon at the top of the Libraries section.

- Navigate to the location of `libSBOLj-2.4.0-withDependencies.jar` on your computer.

- Select the JAR file and click **OK**.

## Usage

To use UCF2SBOL, you’ll need to configure command-line arguments in IntelliJ to specify required parameters. The program accepts multiple arguments, which it uses for SynBioHub login credentials, file paths, and collection details.

### Command-Line Arguments

1. **Open Run/Debug Configurations**:
- Go to **Run** > **Edit Configurations** in IntelliJ.
2. **Set Program Arguments**:
- In the **Program arguments** field, enter the arguments as a single line:

```plaintext
<login email> <password> <login user> <temporary directory> <database prefix> <path to UCF file> [path to UCF input file] [path to UCF output file] [database URL] [collection id] [collection version] [collection name] [collection description] [collection pubMedId]
```

- Example of full commandline arguments for Cello UCF v2 conversion:

```plaintext
"user@example.com" "password" "username" "/tmp" "http://synbiohub.org/" "path/to/ucf.json" "path/to/ucf_input.json" "path/to/ucf_output.json" "http://synbiohub.org/" "Eco1C1G1T1_Parts" "1" "Cello Parts" "These are the Cello parts" "27034378"
```

- Example of the minimal set of commandline arguments for Cello UCF v2 conversion:

```plaintext
"user@example.com" "password" "username" "/tmp" "http://synbiohub.org/" "path/to/ucf.json" "path/to/ucf_input.json" "path/to/ucf_output.json"
```

- Example of the minimal set of commandline arguments for Cello UCF v1 conversion:

```plaintext
"user@example.com" "password" "username" "/tmp" "http://synbiohub.org/" "path/to/ucf.json"
```

- Click **Apply** and **OK** to save your configuration.
3. **Run the Program**:
- To run the program with a sample configuration: Click the **Run** button in IntelliJ to execute the program with your specified command-line arguments.

### Command-Line Arguments Overview

1. **login email** - Email for SynBioHub login.

2. **password** - Password for SynBioHub login.

3. **login user** - Username for SynBioHub login.

4. **temporary directory** - Directory path for storing temporary files.

5. **database prefix** - Prefix for the database URI.

6. **path to UCF file** - Path to the UCF JSON file.

7. **path to UCF input file** (optional) - Path to the input JSON file.

8. **path to UCF output file** (optional) - Path for the output JSON file.

9. **database URL** (optional) - URL for the SynBioHub database. Defaults to the database prefix if not provided.

10. **collection id** (optional) - ID of the collection to use. Defaults to "Eco1C1G1T1_Parts".

11. **collection version** (optional) - Version of the collection. Defaults to "1".

12. **collection name** (optional) - Name of the collection. Defaults to "Cello Parts".

13. **collection description** (optional) - Description of the collection. Defaults to "These are the Cello parts".

14. **collection pubMedId** (optional) - PubMed ID associated with the collection. Defaults to "27034378".
