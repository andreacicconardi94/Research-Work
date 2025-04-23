import pandas as pd
import plotly.express as px
import plotly.graph_objs as go
from dash import Dash, dcc, html, Input, Output
import numpy as np

# Load the dataset
df = pd.read_csv("esempio_ipermappa.txt", delimiter="\t", comment="#", header=None)
df.columns = ["X", "Y", "Raman_Shift", "Intensity"]


# Function to integrate the spectrum at a given (X, Y) point
def integrate_spectrum(group):
    return np.trapezoid(group["Intensity"], group["Raman_Shift"])


# Compute the integral of the spectrum for each spatial point (X, Y)
df_integrated = df.groupby(["X", "Y"]).apply(integrate_spectrum).reset_index(name="Intensity")

# Create Dash app
app = Dash(__name__)
server = app.server  # useful if deploying online later

# App layout
app.layout = html.Div([
    html.H2("Hyperspectral PL Map"),
    dcc.Graph(id='heatmap', config={'displayModeBar': False}),
    html.Div(id='selected-point'),
    dcc.Graph(id='spectrum', style={'width': '750px', 'height': '450px', 'display': 'inline-block'})
])

# Callback to update graphs upon user interaction
@app.callback(
    Output('heatmap', 'figure'),
    Output('selected-point', 'children'),
    Output('spectrum', 'figure'),
    Input('heatmap', 'clickData')
    )

def update_output(clickData):
    # Main PL map using spectrum integration
    heat_fig = px.scatter(df_integrated, x="X", y="Y", color="Intensity",
                          color_continuous_scale="Viridis",
                          symbol_sequence=["square"],
                          title="PL Map (Integrated Intensity)",
                          width=500, height=400)

    heat_fig.update_traces(marker=dict(size=20, line=dict(width=0)))
    heat_fig.update_yaxes(scaleanchor="x", scaleratio=1)

    if clickData is None:
        return heat_fig, "Click a point on the map", go.Figure()

    # Extract clicked coordinates
    clicked_x = clickData["points"][0]["x"]
    clicked_y = clickData["points"][0]["y"]

    # DEBUG PRINT (optional)
    print(f"Clicked: X={clicked_x}, Y={clicked_y}")
    print("First real coordinates in the dataset:")
    print(df[["X", "Y"]].drop_duplicates().head())

    # Filter for the corresponding spectrum
    mask = (df["X"] == clicked_x) & (df["Y"] == clicked_y)
    df_point = df[mask]

    if df_point.empty:
        spec_fig = go.Figure()
        spec_fig.update_layout(
            title="⚠️ No spectrum available for this point",
            xaxis_title="Energy (eV)",
            yaxis_title="Intensity (a.u.)"
        )
        return heat_fig, f"No data for X={clicked_x:.3f}, Y={clicked_y:.3f}", spec_fig

    # Plot the selected spectrum
    spec_fig = go.Figure()
    spec_fig.add_trace(go.Scatter(x=df_point["Raman_Shift"], y=df_point["Intensity"],
                                  mode='lines', line=dict(color='blue')))
    spec_fig.update_layout(
        title=f"Spectrum at X={clicked_x:.3f}, Y={clicked_y:.3f}",
        xaxis_title="Energy (eV)",
        yaxis_title="Intensity (a.u.)"
    )

    return heat_fig, f"Selected: X={clicked_x:.3f}, Y={clicked_y:.3f}", spec_fig

# Run the server
if __name__ == '__main__':
    app.run(debug=True)
