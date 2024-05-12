import pathlib

from ortools.constraint_solver import pywrapcp, routing_enums_pb2

from src.reader import read_data


def create_data_model(matrix_numpy):
    """Stores the data for the problem."""
    data = {"distance_matrix": ([[int(j) for j in i] for i in matrix_numpy]), "num_vehicles": 1, "depot": 0}
    return data


data_numpy, solution_best = read_data(str(pathlib.Path(__file__).parent.parent / 'datasets/data7.txt'))

data = create_data_model(data_numpy)

manager = pywrapcp.RoutingIndexManager(
    len(data["distance_matrix"]), data["num_vehicles"], data["depot"]
)
routing = pywrapcp.RoutingModel(manager)


def distance_callback(from_index, to_index):
    """Returns the distance between the two nodes."""
    # Convert from routing variable Index to distance matrix NodeIndex.
    from_node = manager.IndexToNode(from_index)
    to_node = manager.IndexToNode(to_index)
    return data["distance_matrix"][from_node][to_node]


transit_callback_index = routing.RegisterTransitCallback(distance_callback)
routing.SetArcCostEvaluatorOfAllVehicles(transit_callback_index)

search_parameters = pywrapcp.DefaultRoutingSearchParameters()
search_parameters.first_solution_strategy = (
    routing_enums_pb2.FirstSolutionStrategy.PATH_CHEAPEST_ARC
)
search_parameters.local_search_metaheuristic = (
    routing_enums_pb2.LocalSearchMetaheuristic.GUIDED_LOCAL_SEARCH
)
search_parameters.time_limit.seconds = 30
search_parameters.log_search = True


def print_solution(manager, routing, solution):
    """Prints solution on console."""
    print(f"Objective: {solution.ObjectiveValue()} miles")
    index = routing.Start(0)
    plan_output = "Route for vehicle 0:\n"
    route_distance = 0
    while not routing.IsEnd(index):
        plan_output += f" {manager.IndexToNode(index)} ->"
        previous_index = index
        index = solution.Value(routing.NextVar(index))
        route_distance += routing.GetArcCostForVehicle(previous_index, index, 0)
    plan_output += f" {manager.IndexToNode(index)}\n"
    plan_output += f"Best solution: {solution_best}, Route distance: {route_distance}, {route_distance / solution_best:0.3%}"
    print(plan_output)


solution = routing.SolveWithParameters(search_parameters)
if solution:
    print_solution(manager, routing, solution)
