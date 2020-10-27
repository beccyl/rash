import dash_core_components as dcc
import dash_html_components as html
from dash.dependencies import Input, Output

from app import app
from apps import breast_sample, thyroid_sample, breast_all, thyroid_all
#from apps import breast_sample, breast_all

server = app.server

app.layout = html.Div([
    # represents the URL bar, doesn't render anything
    dcc.Location(id='url', refresh=False),
    html.Div(id='page-content')
])

default_layout = html.Div([
    dcc.Link('Navigate to "/"', href='/'),
    html.Br(),
    dcc.Link('Navigate to "/apps/breast_sample"', href='/apps/breast_sample'),
    html.Br(),
    dcc.Link('Navigate to "/apps/breast_all"', href='/apps/breast_all'),
    html.Br(),
    dcc.Link('Navigate to "/apps/thyroid_sample"', href='/apps/thyroid_sample'),
    html.Br(),
    dcc.Link('Navigate to "/apps/thyroid_all"', href='/apps/thyroid_all'),
])


@app.callback(Output('page-content', 'children'),
              [Input('url', 'pathname')])
def display_page(pathname):
    if pathname == '/apps/breast_sample':
        return breast_sample.layout
    elif pathname == '/apps/breast_all':
        return breast_all.layout
    elif pathname == '/apps/thyroid_sample':
        return thyroid_sample.layout
    elif pathname == '/apps/thyroid_all':
        return thyroid_all.layout
    else:
        return default_layout

        #return '404'

server = app.server

if __name__ == '__main__':
    app.run_server(debug=True)
