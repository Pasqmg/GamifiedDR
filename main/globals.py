"""
Global variable definition.
Adjust the values of these variables to run different experiments and define demand time restrictions.
"""
import os

# Input/Output global variables
INPUT_PATH = "../data/input/"
OUTPUT_PATH = "../data/output/"

# Path to folder containing problem instance files
CONFIG_PATH = "../data/instances/"
# Adjust to routes file name
ROUTES_FILE = os.path.join(INPUT_PATH, 'precomputed_routes.json')
# Adjust to stops file name
STOPS_FILE = os.path.join(INPUT_PATH, 'filtered_stops.json')

# Demand-generation global variables, which affect the time window computation of each Stop within a Request
# OSRM petition url
ROUTE_HOST = "http://localhost:5000/"
# Time-related globals
MAXIMUM_WAITING_TIME_MINUTES = 20
SERVICE_MINUTES_PER_PASSENGER = 1
TRAVEL_FACTOR = 2.5  # to compute maximum on-board time
SLACK_TIME_MINUTES = 5  # minutes of margin for the system to have flexibility
