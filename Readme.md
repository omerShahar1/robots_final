
# description
This project focuses on developing a sophisticated navigation system that utilizes raw Global Navigation Satellite System (GNSS) measurements. The primary objectives include:

Real-time Position Editing: Implementing an efficient and accurate algebraic approach to dynamically adjust and refine real-time positional data.

Satellite Filtering: Incorporating the ability to filter satellites based on several criteria, including:

Constellation: Filtering satellites by their specific GNSS constellation (e.g., GPS, GLONASS, Galileo, BeiDou).
Otzana: Identifying and filtering satellites based on Otzana, a metric related to satellite health and reliability.
False Satellites: Detecting and excluding erroneous satellite signals to ensure data accuracy.
Disruption Handling: Managing and mitigating disruptions, particularly those identified as "Cairo + Beirut" scenarios, which refer to specific types of interference or anomalies that can affect GNSS signals.

Disturbance Identification and Management: Implementing an algorithm to detect disturbances in GNSS signals and effectively address them to maintain the reliability of the navigation data.
