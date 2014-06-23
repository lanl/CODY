function random_pair(upper, ...)
    --[[--

    Returns a pair of random points in the interval [1,upper]. Each
    point is picked with uniform probability. The pair (v0,v1) also has
    the property that v0 != v1. The results are returned in order.  If
    a second argument is passed to the function, it becomes a fixed
    member of the pair.

    --]]--

    local args = table.pack(...)
    local v0, v1

    if args ~= nil and args.n == 1 then
        v0 = args[1]
    else
        v0 = math.random(upper)
    end

    v1 = math.random(upper)

    while v1 == v0 do -- no equality
        v1 = math.random(upper)
    end

    if v0 <= v1 then
        return v0, v1
    else
        return v1, v0
    end

end

function create_graph(graph_type, npoints, nedges, fname)
    --[[--

    This function creates different types of graphs with
    npoints points and
        
        nedges edges for pure random graphs
        nedges * npoints edges for regular random graphs
        npoints * 2 edges for contiguous graphs
        nedges edges for graphs from file

    The graph type must be one of 
        
        'pure_random'
        'regular_random'
        'contiguous'
        'file'

    --]]--

    local t = {
        pure_random    = pure_random_graph,
        regular_random = regular_random_graph,
        contiguous     = regular_contiguous_graph,
        file           = graph_from_file
    }

    if t[graph_type] == nil then
        print("Unknown graph type specified.")
        return {}
    end

    return t[graph_type](npoints, nedges, fname)

end

function pure_random_graph(npoints, nedges) 
    --[[--

    Returns a pure, undirected random graph, i.e. each edge is picked,
    at random, from the interval [1, npoints] with uniform
    probability. We ensure that there are no duplicate edges of the
    form (v0, v1) and (v1, v0) by ordering the vertices in each edge so
    that v0 <= v1.

    --]]--

    -- seed random number generator
    math.randomseed(os.time())

    local e = {}

    for i = 1,nedges do
        v0, v1 = random_pair(npoints)

        while e[v0] ~= nil and e[v0][v1] do -- no redunant edges
            v0, v1 = random_pair(npoints)
        end

        if e[v0] == nil then
            e[v0] = {}
        end

        e[v0][v1] = true
    end

    return e

end

function regular_random_graph(npoints, degree)
    --[[--

    Returns a random semi-regular graph, where each vertex has 
    degree close to 'degree'. Edges are sampled uniformly from 
    the interval 

    [1, npoints] \ v 

    for each vertex v. Note that there are no duplicate edges,
    but the degree of each vertex varies about 'degree'.

    --]]--

    -- seed random number generator
    math.randomseed(17)

    local e = {}

    for i = 1,npoints do
        for j = 1,degree do
            local v0, v1 = random_pair(npoints, i)

            if e[v0] == nil then
                e[v0] = {}
            end

            e[v0][v1] = true            
        end
    end

    return e
end

function regular_contiguous_graph(npoints)
    --[[--

    Returns an undirected graph with two edges per vertex. The
    edges should be contiguous in memory, i.e. edges are of the form
    (i, i+1), (i+1, i+2), (i+2, i+3) etc.

    --]]--

    local e = {}

    for i = 1,npoints do

        local v
        if i < npoints then
            v = i + 1
        else
            v = 1
        end

        if e[i] == nil then
            e[i] = {}
        end

        e[i][v] = true

    end

    return e

end

function graph_to_file(edges, fname)
    --[[--

    Writes the edges in the table 'edges' to the file 'fname', one edge
    per line in the form x,y.

    --]]--

    local file = io.open(fname, "w")

    for k,v in pairs(edges) do 
        for x,_ in pairs(v) do
            file:write(k .. "," .. x .. "\n")
        end
    end

    file:close()

end

function graph_from_file(npoints, nedges, fname)
    --[[--

    Reads and returns a table of edges from the file 'fname'. Edges,
    one per line, should be of the form x,y.

    --]]--

    local e = {}
    local p = {}
    local file = io.open(fname, "r")
    local edge_ct = 0

    for line in file:lines() do
        local v0, v1
        string.gsub(line, "(%d+),(%d+)", function(x,y) v0 = x v1 = y end)

        if e[v0] == nil then
            e[v0] = {}
        end

        e[v0][v1] = true
        p[v0] = true
        p[v1] = true
        edge_ct = edge_ct + 1
    end

    file:close()

    --local point_ct = 0
    --for _ in pairs(p) do point_ct = point_ct + 1 end

    --if edge_ct ~= nedges or point_ct ~= npoints then
    if edge_ct ~= nedges then
        print("Error, graph from file", fname, "does not")
        print("match nedges:", nedges, edge_ct)
        return {}
    end

    return e
end

function test ()

    t0 = os.time()
    --edges = pure_random_graph(5,10)
    --edges = regular_random_graph(10,2)
    --edges = regular_contiguous_graph(10)
    t1 = os.time()
    print("Graph Creation Time: " .. t1 - t0)

    e = 0
    t0 = os.time()
    for k,v in pairs(edges) do
        for x,_ in pairs(v) do
            print(k .. "," .. x)
            e = e + 1
        end
    end
    t1 = os.time()
    print("Graph iteration time: " .. t1 - t0)
    print("Done")
    print("Created graph with " .. e .. " edges")

end



