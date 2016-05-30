
def projection_graph(G, proj_fn, filename=None, verbose=False):
    r"""
    EXAMPLES::

        sage: from slabbe.bispecial_extension_type import rec_enum_set_under_language_joined_from_pairs
        sage: R = rec_enum_set_under_language_joined_from_pairs(pairs, LBrun, S, keep_empty=False, label='previous', growth_limit=1)
        sage: from slabbe.bispecial_extension_type import recursively_enumerated_set_to_digraph
        sage: %time G = recursively_enumerated_set_to_digraph(R)
        CPU times: user 15min 59s, sys: 5.07 s, total: 16min 4s
        Wall time: 35min 17s
        sage: G
        Looped multi-digraph on 8341 vertices

    ::

        sage: f2 = lambda ext:(ext.left_valence(),len(ext.left_word_extensions()))
        sage: f3 = lambda ext:(ext.left_valence(),len(ext.left_word_extensions()),ext.multiplicity())
        sage: f4 = lambda ext:(ext.left_valence(),len(ext.left_word_extensions()),ext.multiplicity(),ext.is_ordinaire())
        sage: map_f2 = lambda s:tuple(sorted(map(f2,s)))
        sage: map_f3 = lambda s:tuple(sorted(map(f3,s)))
        sage: map_f4 = lambda s:tuple(sorted(map(f4,s)))
        sage: projection_graph(G, map_f2, key='2')

    ::

        sage: from slabbe.mult_cont_frac import Brun
        sage: algo = Brun()
        sage: S = algo.substitutions()
        sage: from slabbe.language import languages
        sage: LBrun = languages.Brun()
        sage: E1, E2, E3, E4, E5 = empty_word_extension_types()
        sage: pairs = [(E1, 312), (E2, 312), (E3, 312), (E4, 321), (E5, 321)]
        sage: f = lambda S:any(len(ext.left_word_extensions())>2 for ext in S)
        sage: from slabbe.bispecial_extension_type import rec_enum_set_under_language_joined_from_pairs
        sage: R = rec_enum_set_under_language_joined_from_pairs(pairs, LBrun, S, keep_empty=False, label='previous', growth_limit=1, filter_fn=f)
        sage: R
        A recursively enumerated set (breadth first search)
        sage: from slabbe.bispecial_extension_type import recursively_enumerated_set_to_digraph
        sage: G = recursively_enumerated_set_to_digraph(R)
        sage: G
        Looped multi-digraph on 453 vertices
        sage: 2216+2541+1796+1220+120+453
        8346
        sage: projection_graph(G, map_f4, '4')
         ((3, 5, 0, False),): 3, 
         ((3, 5, 0, True),): 2, 
         ((2, 4, 0, True),): 3,
         ((3, 4, 0, False),): 2, 
         ((2, 3, 0, True), (2, 3, 0, True)): 8}

         ((2, 3, 0, True),): 14, 
         ((3, 3, 0, False),): 18, 
         ((3, 3, 0, True),): 12, 
        {((2, 2, 0, True), (2, 3, 0, True)): 157, 
         ((2, 2, 0, True), (3, 3, 0, True)): 12, 
         ((2, 2, 0, True), (3, 3, 0, False)): 12, 
         ((2, 2, -1, False), (2, 3, 1, False)): 43, 
         ((2, 2, 0, True), (2, 2, 0, True), (2, 3, 0, True)): 131, 
         ((2, 2, -1, False), (2, 2, 0, True), (2, 3, 1, False)): 36, 
    """
    edges = set((proj_fn(A[0]),proj_fn(B[0])) for A,B,_ in G.edges())
    G_proj = DiGraph(edges, format='list_of_edges', loops=True, multiedges=True)
    if verbose:
        d = dict(Counter(proj_fn(s[0]) for s in G.vertices()))
        rows = [(value, key) for key,value in d.iteritems()]
        rows.sort(reverse=True,key=lambda row:row[1])
        print table(rows=rows)
    if filename:
        print TikzPicture.from_graph(G_proj, prog='dot').pdf(filename)
    return G_proj

def digraph_move_label_to_edge(G, label_function=None):
    r"""
    """
    if label_function:
        edges = [(u,v,label_function(label)) for ((u,_), (v,label), _) in G.edges()]
    else:
        edges = [(u,v,label) for ((u,_), (v,label), _) in G.edges()]
    return DiGraph(edges, format='list_of_edges', loops=True, multiedges=True)

def induced_subgraph(G, filter):
    r"""
    Removes vertices by keeping the adjacenjcies.

    TODO:

        simplify the edges aaaaa*aaa*aaa* to a* only
    """
    GG = G.copy()
    loops = dict((u, label) for (u,v,label) in GG.loop_edges())
    for v in GG.vertices():
        if filter(v):
            continue
        incoming = [(x,y,l) for (x,y,l) in GG.incoming_edges(v) if x != y]
        outgoing = [(x,y,l) for (x,y,l) in GG.outgoing_edges(v) if x != y]
        GG.delete_vertex(v)
        if v in loops:
            it = itertools.product(outgoing, [(loops[v][0]*10,)], incoming) 
        else:
            it = itertools.product(outgoing, [tuple()], incoming) 
        for c,b,a in it:
            _,y,labelout = c
            x,_,labelin = a
            GG.add_edge(x,y, labelout + b + labelin)
    return GG

