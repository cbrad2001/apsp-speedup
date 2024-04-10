import sys
import random

def generate_connected_graph(num_nodes, num_edges):
    if num_edges < num_nodes - 1:
        print("Error: Number of edges must be at least num_nodes - 1 to create a connected graph.")
        return
    if num_edges > (num_nodes * (num_nodes - 1)):
        print("Error: Too many edges.")
        return

    # Initialize the graph with each node connected to at least one other node
    edges = [(i, (i + 1) % num_nodes) for i in range(num_nodes)]

    # Generate additional random edges
    for _ in range(num_edges - (num_nodes)):
        source = random.randint(0, num_nodes - 1)
        target = random.randint(0, num_nodes - 1)
        while target == source or (source, target) in edges:
            target = random.randint(0, num_nodes - 1)
        edges.append((source, target))

    # Write edges to a text file
    with open("testGraph.txt", "w") as file:
        for edge in edges:
            file.write("{} {}\n".format(edge[0], edge[1]))

if __name__ == "__main__":
    if len(sys.argv) != 3:
        print("Usage: python generate_graph.py <num_nodes> <num_edges>")
        sys.exit(1)

    num_nodes = int(sys.argv[1])
    num_edges = int(sys.argv[2])

    generate_connected_graph(num_nodes, num_edges)
