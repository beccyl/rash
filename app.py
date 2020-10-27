import dash
import flask

#server = flask.Flask(__name__) # define flask app.server
#app = dash.Dash(__name__, server=server, suppress_callback_exceptions = True) # call flask server
#app.config.suppress_callback_exceptions = True
app = dash.Dash(__name__, suppress_callback_exceptions = True) # call flask server
#server = app.server
