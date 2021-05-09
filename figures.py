"""
Caleb Ellington, Anush Devadhasan, Daniel Jeong
"""

from IPython.display import Image
from causalnex.plots import plot_structure, NODE_STYLE, EDGE_STYLE


def visualize_dbn(dbn, savepath=None):
    viz = plot_structure(
        dbn,
        graph_attributes={"scale": "2.0"},
        all_node_attributes=NODE_STYLE.NORMAL,
        all_edge_attributes=EDGE_STYLE.NORMAL)
    img = viz.draw(format='png')
    if savepath is not None:
        with open(savepath, 'wb') as file:
            file.write(img)
    return img
#     Image(viz.draw(format='png'))
