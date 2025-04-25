import heapq
import math

import numpy
import numpy as np

from main import database
from main.insertion import Insertion
from main.itinerary import Itinerary
from main.stop import Stop


class Scheduler:
    """
    Scheduler object. The Scheduler creates and solves a Demand-responsive problem instance,
    created according to the input data contained in the Database.

    As of the current version, the Scheduler launches a one-shot process that iteratively schedules
    customer requests to vehicle itineraries. Request scheduling can be performed by:
        - issuance time (see schedule_all_requests_by_time_order)
        - minimal cost (see schedule_all_requests_by_minimal_cost)

    Requests define a customer trip. The assignment of a request's trip to a vehicle itinerary implies
    the preservation of the trip's spatial and temporal constraints. Requests are assigned to the best
    possible itinerary following the objective function defined in itinerary.compute_cost(). The search
    for the best itinerary is done through exhaustive search of all feasible itinerary positions where
    a given customer trip can be inserted. Such best position defined the 'best insertion', the object
    matching a request with its assigned itinerary and position within it.

    A Demand-responsive problem instance is solved once all requests have been either scheduled to a
    vehicle itinerary or rejected.
    """

    def __init__(self, database):
        """
        The Scheduler initialises the Database object of the given problem instance.
        The Scheduler manages
            # lists of pending, scheduled and rejected requests
            # vehicle itineraries
            # problem solution metrics (simulation_dict)
        """
        self.db = database
        # List of unscheduled Requests, ordered by time of issuance
        self.pending_requests = []
        # List of accepted Requests
        self.scheduled_requests = []
        # List of rejected Requests
        self.rejected_requests = []
        # List of fleet vehicles
        self.vehicles = []
        # List of itineraries (initially empty)
        self.itineraries = []
        # Dictionary of insertions of each itinerary (initially empty)
        self.itinerary_insertion_dic = {}
        # Dictionary wit the evaluation metrics of a problem solution
        self.simulation_dict = {}

    ################################################
    ############# Auxiliary methods ################
    ################################################
    def delete_pending_request(self, passenger_id):
        self.pending_requests = [x for x in self.pending_requests if x.passenger_id != passenger_id]

    ################################################
    ########## Insertion search methods ############
    ################################################

    def exhaustive_search(self, t, verbose=0):
        # list to store the found insertions
        feasible_insertions = []
        # Minimum cost increment found for the insertion of the current request
        min_delta = math.inf
        # Insertion with the minimum cost increment found so far
        best_insertion = None

        # Extract Request's stops
        Spu = new_stop_from_stop(t.Spu)
        Ssd = new_stop_from_stop(t.Ssd)
        passenger_id = t.passenger_id

        if verbose > 0:
            print("Searching for best insertion for trip {}".format(t.passenger_id))
        # for each vehicle
        for I in self.itineraries:
            if verbose > 0:
                print("\tSearching inside itinerary {}".format(I.vehicle_id))
            # Filter list of stops to keep only those not yet visited
            index_current = I.stop_list.index(I.current_loc)
            filtered_stops_i = [new_stop_from_stop(x) for x in I.stop_list[index_current:]]
            # Find feasible insertion for Spu
            for index_stop_i in range(len(filtered_stops_i) - 1):
                if verbose > 0:
                    print("\t\tTesting insertion of Spu in position {}".format(index_stop_i + index_current + 1))
                # extract leg R -> T
                # DEBUG
                try:
                    R = filtered_stops_i[index_stop_i]
                except IndexError:
                    print("ERROR Searching inside itinerary {}".format(I.vehicle_id))
                    print(I.to_string())
                    print()
                    print("with the following list of filtered stops: {}".format([x.id for x in filtered_stops_i]))
                    for x in filtered_stops_i:
                        print(x.to_string())
                    print()
                    print("and an index_stop_i of: {}".format(index_stop_i))
                    print(len(filtered_stops_i), index_stop_i)
                    exit()

                T = R.snext
                # Check feasibility of inserting Spu in R's position, so that leg (R -> R.rnext)
                # becomes (Spu -> R.snext) therefore creating also a new leg (R -> Spu)
                test, code = I.pickup_insertion_feasibility_check(t, Spu, R, T)
                if test:
                    if verbose > 0:
                        print("\t\t\tfeasible")
                    # Once we select a feasible leg to insert Spu, store the index
                    index_Spu = index_stop_i + index_current + 1
                    # Copy of the itinerary to avoid modifications over the original
                    I_with_Spu = new_itinerary_from_itinerary(I)
                    # I_with_Spu = copy_Itinerary(I)
                    # Insert Spu in the itinerary and re-calculate EAT carried forward over its putative successors
                    I_with_Spu.insert_stop(Spu, index_Spu)
                    # Compute the insertion's net additional cost
                    I_with_Spu.compute_cost()
                    delta_i = I_with_Spu.cost - I.cost
                    # If net additional cost < minimum cost increment found so far, go on to insert Ssd
                    if delta_i < min_delta:
                        # Look for a leg to insert Ssd in each stop in the itinerary after R

                        # Filter list of stops to keep only those not yet visited
                        filtered_stops_j = [new_stop_from_stop(x) for x in I_with_Spu.stop_list[index_Spu:]]

                        for index_stop_j in range(len(filtered_stops_j) - 1):
                            if verbose > 0:
                                print("\t\t\t\tTesting insertion of Ssd in position {}"
                                      .format(index_stop_j + index_Spu + 1))
                            R = filtered_stops_j[index_stop_j]
                            T = R.snext
                            test, code = I_with_Spu.setdown_insertion_feasibility_check(t, index_Spu,
                                                                                        index_stop_j + index_Spu + 1,
                                                                                        I_with_Spu.stop_list, Ssd, R, T)
                            if test:
                                if verbose > 0:
                                    print("\t\t\t\t\tfeasible")
                                # Once we select a feasible leg to insert Ssd, store the index
                                index_Ssd = index_stop_j + index_Spu + 1
                                # Copy of the itinerary to avoid modifications over the original
                                I_with_Spu_Ssd = new_itinerary_from_itinerary(I_with_Spu)
                                # I_with_Spu_Ssd = copy_Itinerary(I_with_Spu)
                                # I_with_Spu_Ssd = copy_Itinerary(I_with_Spu)
                                I_with_Spu_Ssd.insert_stop(Ssd, index_Ssd)
                                # Compute the insertion's net additional cost
                                I_with_Spu_Ssd.compute_cost()
                                delta_ij = I_with_Spu_Ssd.cost - I.cost

                                # Create insertion object and store it in the list
                                found = Insertion(itinerary=I, trip=t, index_Spu=index_Spu, index_Ssd=index_Ssd,
                                                  cost_increment=delta_ij)
                                feasible_insertions.append((found, delta_ij))

                                # if delta_ij < minimum cost increment found so far, update minimum cost
                                if delta_ij < min_delta:
                                    min_delta = delta_ij
                                    best_insertion = found
                            else:
                                if verbose > 0:
                                    print("\t\t\t\tunfeasible")
                                # Try in next Spu's position
                                if code == 0:
                                    break
                        # end of filtered_stops_j for
                    # end of delta_i < delta_min check
                else:
                    if verbose > 0:
                        print("\t\tunfeasible")
                    # Go to next itinerary
                    if code == 0:
                        break
            # end of filtered_stops_i for

        if verbose > 0:
            print()
        return best_insertion, feasible_insertions

    def get_minimal_cost_insertion(self, verbose=0):
        found_insertions = []
        for request in self.pending_requests:
            best_insertion, feasible_insertions = self.exhaustive_search(request)
            if verbose > 0:
                print("Found {} feasible insertion(s)".format(len(feasible_insertions)))
            for x in feasible_insertions:
                found_insertions.append(x)
        if len(found_insertions) > 0:
            found_insertions.sort(key=lambda a: a[1])
            if verbose > 0:
                print("Found insertions: {}".format(found_insertions))
            return found_insertions[0]
        else:
            return None, None

    ################################################
    ########### Trip operation methods #############
    ################################################

    def insert_trip(self, insertion):

        # Add insertion to the itinerary_insertion_dic of the Scheduler
        vehicle_id = insertion.I.vehicle_id
        passenger_id = insertion.t.passenger_id
        self.itinerary_insertion_dic[vehicle_id].append(insertion)

        # Extract Request attributes
        Spu, Ssd = insertion.t.Spu, insertion.t.Ssd
        Spu.passenger_id = passenger_id
        Ssd.passenger_id = passenger_id

        # Insert pickup stop and update itinerary times
        insertion.I.insert_stop(Spu, insertion.index_Spu)

        # Insert setdown stop and update itinerary times
        insertion.I.insert_stop(Ssd, insertion.index_Ssd)

        # Initial definition of passenger-loading variables
        # S = Spu
        # S.npass = S.sprev.npass
        insertion.I.stop_list[insertion.index_Spu].npass = insertion.I.stop_list[insertion.index_Spu - 1].npass
        # S.npres = S.sprev.npres
        insertion.I.stop_list[insertion.index_Spu].npres = insertion.I.stop_list[insertion.index_Spu - 1].npres
        # S = Ssd
        # S.npass = S.sprev.npass
        insertion.I.stop_list[insertion.index_Ssd].npass = insertion.I.stop_list[insertion.index_Ssd - 1].npass
        # S.npres = S.sprev.npres
        insertion.I.stop_list[insertion.index_Ssd].npres = insertion.I.stop_list[insertion.index_Ssd - 1].npres

        # Adjust passenger-loading variables for stops between Spu and Ssd.sprev (both included)
        npshare_t = insertion.I.capacity - insertion.t.npass
        for i in range(insertion.index_Spu, insertion.index_Ssd):
            insertion.I.stop_list[i].npass = insertion.I.stop_list[i].npass + insertion.t.npass
            insertion.I.stop_list[i].npres = insertion.I.stop_list[i].npres + (insertion.I.capacity - npshare_t)

        # Update itinerary distance and time cost
        insertion.I.traveled_km = insertion.I.compute_traveled_km()
        insertion.I.cost = insertion.I.compute_cost()

        # check for duplicate stops
        # insertion.I.check_duplicated_stops()

    def remove_trip(self, insertion, verbose=0):
        # Extract Request attributes
        Spu, Ssd = insertion.t.Spu, insertion.t.Ssd

        if verbose > 0:
            print("Removing stops in indexes {} and {}\n".format(insertion.index_Spu, insertion.index_Ssd))
            print("\tthese are stops {} and {}\n"
                  .format(insertion.I.stop_list[insertion.index_Spu].id, insertion.I.stop_list[insertion.index_Ssd].id))
            print("The stops between the removed trip are {}\n"
                  .format([x.id for x in insertion.I.stop_list[insertion.index_Spu + 1:insertion.index_Ssd]]))

        # Remove setdown stop and update itinerary times
        insertion.I.remove_stop(Ssd, insertion.index_Ssd)

        # Remove pickup stop and update itinerary times
        insertion.I.remove_stop(Spu, insertion.index_Spu)

        if insertion.index_Ssd - insertion.index_Spu > 1:
            # New index for Spu.snext is index_Spu
            i = insertion.index_Spu
            # New index for Ssd.sprev is (index_Ssd - index_Spu)
            j = insertion.index_Ssd - insertion.index_Spu + 1

            if verbose > 0:
                print("Now that the trip has been removed, the stops that were between Spu and Ssd have new indexes..."
                      "\n\tthese are comprised between {} and {}, both included\n".format(i, j))
                print("Let's check it:\n\t Stop in index {}: {}\n\tStop in index {}: {}"
                      .format(i, insertion.I.stop_list[i].id, j, insertion.I.stop_list[j].id))
                print("The whole list between {} and {} is: {}\n"
                      .format(i, j, [x.id for x in insertion.I.stop_list[i:j + 1]]))
                print("The final stop list is {}\n".format([x.id for x in insertion.I.stop_list]))

            npshare_t = insertion.I.capacity - insertion.t.npass
            for index_S in range(i, j + 1):
                insertion.I.stop_list[index_S].npass = insertion.I.stop_list[index_S].npass - insertion.t.npass
                insertion.I.stop_list[index_S].npres = insertion.I.stop_list[index_S].npres - (
                        insertion.I.capacity - npshare_t)

        # Update itinerary distance and time cost
        insertion.I.traveled_km = insertion.I.compute_traveled_km()
        insertion.I.cost = insertion.I.compute_cost()

    ################################################
    ############# Scheduling methods ###############
    ################################################

    def schedule_all_requests_by_minimal_cost(self, verbose=0):
        max_tries = len(self.pending_requests) * 5
        counter = 0
        while len(self.pending_requests) > 0 and counter < max_tries:
            next_insertion, ci = self.get_minimal_cost_insertion()
            if next_insertion is not None:
                if verbose > 0:
                    print("Scheduling best insertion")
                    print(next_insertion.to_string())
                self.insert_trip(next_insertion)
                self.delete_pending_request(next_insertion.t.passenger_id)
            counter += 1
        if verbose > 0:
            if counter < max_tries:
                print("All requests scheduled.")
            else:
                for req in self.pending_requests:
                    print("Pending request {}".format(req.passenger_id))
                print()
            self.print_itineraries()

    def schedule_all_requests_by_time_order(self, verbose=0):
        pending_req = len(self.pending_requests)
        while len(self.pending_requests) > 0:
            t = self.pending_requests[0]
            self.pending_requests = self.pending_requests[1:]
            best_insertion, _ = self.exhaustive_search(t)
            # best_insertion, _ = self.exhaustive_search_inplace(t)
            if best_insertion is None:
                if verbose > 1:
                    print("Trip {} can not be scheduled".format(t.passenger_id))
                self.rejected_requests.append(t)
            else:
                self.insert_trip(best_insertion)
                self.scheduled_requests.append(t)
                if verbose > 1:
                    print(best_insertion.to_string())
        for I in self.itineraries:
            # I.update_time_windows()
            I.compute_dispatching()
        rejected_req = len(self.rejected_requests)
        scheduled_req = pending_req - rejected_req
        if verbose > 0:
            self.print_itineraries()
            if verbose > 1:
                self.print_itineraries_debug()
            if verbose > 2:
                if len(self.rejected_requests) > 0:
                    for req in self.rejected_requests:
                        print("Rejected request {}".format(req.passenger_id))
            print("Scheduled requests: {:3d}, Accepted rate: {:.2f}%".format(scheduled_req,
                                                                             (scheduled_req / pending_req) * 100))
            print("Rejected requests: {:3d}, {:.2f}%".format(rejected_req, (rejected_req / pending_req) * 100))
            global_customer_waitings = []
            local_avgs = []
            local_stds = []
            for I in self.itineraries:
                for w in I.customer_waitings:
                    global_customer_waitings.append(w)
                local_customer_waitings = numpy.array(I.customer_waitings)

                local_avg = np.average(local_customer_waitings)
                local_avgs.append(local_avg)

                local_std = np.std(local_customer_waitings)
                local_stds.append(local_std)

            global_customer_waitings = numpy.array(global_customer_waitings)
            print("Customer waiting time: Sum of waitings: {:.2f} min, avg: {:.2f} min, std: {:.2f} min\n".format(
                numpy.sum(global_customer_waitings), numpy.average(global_customer_waitings),
                numpy.std(global_customer_waitings)))

            print("{:.2f}\t{:.2f}\t{:.2f}".format((scheduled_req / pending_req) * 100,
                                                  numpy.average(global_customer_waitings),
                                                  numpy.std(global_customer_waitings)))
            print(self.simulation_stats())

    ################################################
    ######### Solution evaluation methods ##########
    ################################################

    def simulation_stats(self):
        """
        Returns a dictionary grouping all metrics that evaluate a problem solution. These include:
            # overall performance metrics (service quality and costs)
            # for each fleet vehicle
                - vehicle's itinerary and performance metrics
                - for each request assigned to the vehicle's itinerary
                    -- request's metrics
        """
        simulation_dict = {}
        # number of requests
        simulation_dict['total_requests'] = len(self.scheduled_requests) + len(self.rejected_requests)
        # scheduled requests and percentage
        simulation_dict['scheduled_requests'] = len(self.scheduled_requests)
        simulation_dict['scheduled_percent'] = len(self.scheduled_requests) / simulation_dict['total_requests'] * 100
        # rejected requests and percentage
        simulation_dict['rejected_requests'] = len(self.rejected_requests)
        simulation_dict['rejected_percent'] = len(self.rejected_requests) / simulation_dict['total_requests'] * 100
        # total cost and kms
        simulation_dict['total_cost'] = -1
        simulation_dict['total_kms'] = -1
        # number of vehicles
        simulation_dict['num_vehicles'] = len(self.itineraries)
        simulation_dict['vehicle_stats'] = {}
        total_cost = 0
        total_kms = 0
        # vehicle_stats indexed by vehicle_id
        for I in self.itineraries:
            simulation_dict['vehicle_stats'][I.vehicle_id] = I.vehicle_stats()
            total_cost += simulation_dict['vehicle_stats'][I.vehicle_id]['cost']
            total_kms += simulation_dict['vehicle_stats'][I.vehicle_id]['itinerary_kms']
            # vehicle_stats dict contains the itinerary's customer_dict
            simulation_dict['vehicle_stats'][I.vehicle_id]['customer_stats'] = I.customer_stats()
        simulation_dict['total_cost'] = total_cost
        simulation_dict['total_kms'] = total_kms

        self.simulation_dict = simulation_dict
        return simulation_dict

    ################################################
    ######## Visualisation & Debug methods #########
    ################################################

    def print_itineraries(self):
        print("\nPrinting all itineraries:\n")
        for I in self.itineraries:
            print(I.to_string() + "\n")
            print(I.customer_stats())
            print(I.vehicle_stats())

    def print_itineraries_debug(self):
        print("\nPrinting all itineraries:\n")
        for I in self.itineraries:
            print(I.to_string_debug() + "\n")

    def print_itineraries_detail(self):
        for I in self.itineraries:
            print(I.to_string())
            for stop in I.stop_list:
                print(stop.to_string())
            print()

    def print_itinerary_insertion_dic(self):
        print("\nPrinting itinerary-insertion dic:\n")
        for key in self.itinerary_insertion_dic.keys():
            insertion_list = self.itinerary_insertion_dic[key]
            print("\tItinerary of vehicle {} contains {:3d} insertions...\n".format(key, len(insertion_list)))
            for insertion in insertion_list:
                print(insertion.to_string())

    def test_insertion(self):
        print("Initial itinerary ---------------------\n")
        itinerary1 = self.itineraries[0]
        for stop in itinerary1.stop_list:
            print(stop.to_string())
        # print(itinerary1.to_string_debug())
        request1 = self.pending_requests[0]
        insertion1 = Insertion(itinerary=itinerary1, trip=request1, index_Spu=1, index_Ssd=2, cost_increment=10)
        print("Itinerary after insertion:-----------------\n")
        self.insert_trip(insertion1)
        for stop in itinerary1.stop_list:
            print(stop.to_string())
        print(itinerary1.to_string_debug())

    # End of Scheduler class


################################################
#### Functions for safe object duplication #####
################################################

def new_stop_from_stop(S):
    # Database
    db = S.db
    # Stop-independent attributes (stop S)
    id = S.id
    # spatial location of S
    coords = S.coords

    # Trip-dependent attributes
    # start of the time-window of S
    start_time = S.start_time
    # end of the time-window of S
    end_time = S.end_time
    # duration of service time at S, normally for loading/unloading passengers
    service_time = S.service_time
    # latest feasible arrival time at S
    lat = S.latest

    # Itinerary-dependent attributes (itinerary I, stop S)
    # predecessor stop in I
    sprev = S.sprev
    # successor stop in I
    snext = S.snext
    # number of passengers on board the vehicle on departure from S
    npass = S.npass
    # number of seats reserved on departure from S
    npres = S.npres
    # CAPACITY CONSTRAINT: capacity <= npass + npress
    # The number of passengers on board + the number of reserved seats can not exceed the capacity of the
    # vehicle represented by the itenerary to which the stop is assigned
    # travel time on the leg connecting S to its successor T
    leg_time = S.leg_time
    # Earliest arrival time
    eat = S.eat
    # Latest departure time
    ldt = S.ldt
    # Earliest feasible commencement of a service interval at S
    eat_f = S.eat_f
    # Latest feasible termination of a service interval at S
    ldt_f = S.ldt_f
    # Slack time
    slack = S.slack

    # Dispatching strategy dependent attributes
    arrival_time = S.arrival_time
    departure_time = S.departure_time

    new_S = Stop(db, id)
    new_S.start_time = start_time
    new_S.end_time = end_time
    new_S.service_time = service_time
    new_S.latest = end_time - service_time
    new_S.sprev = sprev
    new_S.snext = snext
    new_S.npass = npass
    new_S.npres = npres
    new_S.leg_time = leg_time
    new_S.eat = eat
    new_S.ldt = ldt
    new_S.eat_f = eat_f
    new_S.ldt_f = ldt_f
    new_S.slack = slack
    new_S.arrival_time = arrival_time
    new_S.departure_time = departure_time
    return new_S


def new_itinerary_from_itinerary(I):
    db = I.db
    # vehicle to which the itinerary I is assigned
    vehicle_id = I.vehicle_id
    # Capacity of the vehicle
    cap = I.capacity
    # Stop at which the vehicle is located at the beginning of its shift
    start_stop_id = I.start_stop.id
    # Time at which the vehicle begins its shift
    start_time = I.start_time
    # Stop at which the vehicle must be located at the end of its shift
    end_stop_id = I.end_stop.id
    # Time at which the vehicle ends its shift
    end_time = I.end_time
    # List of Stop objects that constitutes the itinerary I
    stop_list = []
    for stop in I.stop_list:
        # stop_list.append(stop)
        stop_list.append(new_stop_from_stop(stop))

    new_I = Itinerary(db, vehicle_id, cap, start_stop_id, end_stop_id,
                      start_time, end_time)
    new_I.stop_list = stop_list
    new_I.compute_traveled_km()
    new_I.compute_cost()

    return new_I
