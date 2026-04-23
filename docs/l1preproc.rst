Level-1 Pre-Processor
=====================


Schematic
---------

.. mermaid::

    graph TD
        CFG("Configuration") --> |"init()"| A["Level-1 PreProcessor"]
        A -->|"process_input_files()"| A1[Clear l1 stack]
        A1 --> B(("Main Loop <br> Source files"))
        B -->|"input_adapter.get_l1()"| C("Read source file")
        C ---|"stage: post_source"| D("Processor items")
        D --> E("Extract polar <br> ocean segments")
        E -->|"stage: post_polar_ocean_segment_extraction" | F("Processor items")
        F --> F1(Add segments to l1 stack)
        F1 -->G{"All l1 segments <br> connected?"}
        G -->|yes| B
        G -->|no| H("Merge connected <br> l1 segments")
        H --> H1["Merged l1 segments"]
        H1 --- |"stage: post_merge"| I("Processor items")
        I --> J("Export l1p netCDF")
        H --> H2["Remaining l1 segments"]
        H2 --> H3["Reset l1 stack"]
        H3 --> B
