from pyvis.network import Network


def vis(G):
    gvis = Network(height='700px', width='700px')
    gvis.from_nx(G)

    gvis.toggle_physics(True)
    gvis.show_buttons()
    gvis.show("Visualization of INFECTTIOUS daily graph.html")

