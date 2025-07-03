import dash
from dash import dcc, html, Input, Output, callback
import plotly.graph_objects as go
import pandas as pd
import yfinance as yf
from datetime import datetime, timedelta
import requests
from bs4 import BeautifulSoup


def get_sp500_symbols():
    """
    Fetch S&P 500 stock symbols from Wikipedia
    """
    try:
        url = 'https://en.wikipedia.org/wiki/List_of_S%26P_500_companies'
        response = requests.get(url)
        soup = BeautifulSoup(response.content, 'html.parser')
        
        # Find the table with S&P 500 companies
        table = soup.find('table', {'class': 'wikitable sortable'})
        symbols = []
        
        for row in table.find_all('tr')[1:]:  # Skip header row
            cells = row.find_all('td')
            if cells:
                symbol = cells[0].text.strip()
                company_name = cells[1].text.strip()
                symbols.append({'label': f"{symbol} - {company_name}", 'value': symbol})
        
        return symbols
    except Exception as e:
        print(f"Error fetching S&P 500 symbols: {e}")
        # Fallback to a subset of popular S&P 500 stocks
        return [
            {'label': 'AAPL - Apple Inc.', 'value': 'AAPL'},
            {'label': 'MSFT - Microsoft Corporation', 'value': 'MSFT'},
            {'label': 'GOOGL - Alphabet Inc.', 'value': 'GOOGL'},
            {'label': 'AMZN - Amazon.com Inc.', 'value': 'AMZN'},
            {'label': 'TSLA - Tesla Inc.', 'value': 'TSLA'},
            {'label': 'META - Meta Platforms Inc.', 'value': 'META'},
            {'label': 'NVDA - NVIDIA Corporation', 'value': 'NVDA'},
            {'label': 'JPM - JPMorgan Chase & Co.', 'value': 'JPM'},
            {'label': 'JNJ - Johnson & Johnson', 'value': 'JNJ'},
            {'label': 'V - Visa Inc.', 'value': 'V'},
        ]


def fetch_stock_data(symbol, period='1y'):
    """
    Fetch stock data for a given symbol from Yahoo Finance
    """
    try:
        stock = yf.Ticker(symbol)
        hist = stock.history(period=period)
        
        if hist.empty:
            return None
        
        # Reset index to make Date a column
        hist.reset_index(inplace=True)
        
        # Create the data structure as requested (date, stock, price)
        data = []
        for _, row in hist.iterrows():
            data.append({
                'date': row['Date'],
                'stock': symbol,
                'price': row['Close']
            })
        
        return pd.DataFrame(data)
    
    except Exception as e:
        print(f"Error fetching data for {symbol}: {e}")
        return None


# Initialize the Dash app
app = dash.Dash(__name__)

# Get S&P 500 symbols
sp500_symbols = get_sp500_symbols()

# Define the app layout
app.layout = html.Div([
    html.H1("S&P 500 Stock Dashboard", 
            style={'textAlign': 'center', 'marginBottom': 30, 'color': '#2c3e50'}),
    
    html.Div([
        html.Label("Select a Stock:", 
                  style={'fontSize': 16, 'fontWeight': 'bold', 'marginBottom': 10}),
        dcc.Dropdown(
            id='stock-dropdown',
            options=sp500_symbols,
            value='AAPL',  # Default selection
            style={'marginBottom': 20}
        ),
    ], style={'width': '50%', 'margin': 'auto', 'marginBottom': 30}),
    
    html.Div([
        dcc.Graph(id='stock-chart')
    ], style={'margin': '20px'}),
    
    html.Div([
        html.H3("Stock Data Table", 
                style={'textAlign': 'center', 'marginTop': 30, 'color': '#2c3e50'}),
        html.Div(id='stock-table')
    ], style={'margin': '20px'})
])


@app.callback(
    [Output('stock-chart', 'figure'),
     Output('stock-table', 'children')],
    [Input('stock-dropdown', 'value')]
)
def update_dashboard(selected_stock):
    """
    Update the chart and table when a stock is selected
    """
    if not selected_stock:
        return {}, "No stock selected"
    
    # Fetch stock data
    df = fetch_stock_data(selected_stock)
    
    if df is None or df.empty:
        return {
            'data': [],
            'layout': {
                'title': f'No data available for {selected_stock}',
                'xaxis': {'title': 'Date'},
                'yaxis': {'title': 'Price ($)'}
            }
        }, "No data available"
    
    # Create the price chart
    fig = go.Figure()
    
    fig.add_trace(go.Scatter(
        x=df['date'],
        y=df['price'],
        mode='lines',
        name=f'{selected_stock} Price',
        line=dict(color='#3498db', width=2)
    ))
    
    fig.update_layout(
        title=f'{selected_stock} Stock Price - Last Year',
        xaxis_title='Date',
        yaxis_title='Price ($)',
        template='plotly_white',
        height=500,
        hovermode='x unified'
    )
    
    # Create a summary table
    latest_price = df['price'].iloc[-1]
    highest_price = df['price'].max()
    lowest_price = df['price'].min()
    price_change = latest_price - df['price'].iloc[0]
    price_change_pct = (price_change / df['price'].iloc[0]) * 100
    
    table_data = html.Table([
        html.Thead([
            html.Tr([
                html.Th('Metric', style={'padding': '10px', 'backgroundColor': '#f8f9fa'}),
                html.Th('Value', style={'padding': '10px', 'backgroundColor': '#f8f9fa'})
            ])
        ]),
        html.Tbody([
            html.Tr([html.Td('Current Price'), html.Td(f'${latest_price:.2f}')]),
            html.Tr([html.Td('52-Week High'), html.Td(f'${highest_price:.2f}')]),
            html.Tr([html.Td('52-Week Low'), html.Td(f'${lowest_price:.2f}')]),
            html.Tr([html.Td('Year-to-Date Change'), html.Td(f'${price_change:.2f}')]),
            html.Tr([html.Td('Year-to-Date % Change'), 
                    html.Td(f'{price_change_pct:.2f}%', 
                           style={'color': 'green' if price_change_pct > 0 else 'red'})]),
        ])
    ], style={'margin': 'auto', 'width': '50%', 'border': '1px solid #ddd'})
    
    return fig, table_data


if __name__ == '__main__':
    app.run_server(debug=True, host='0.0.0.0', port=8050)