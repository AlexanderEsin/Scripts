#!/usr/bin/env python3
"""
Example script demonstrating how to fetch S&P 500 stock data
in the format requested: each row has date, stock, and price.
"""

import pandas as pd
import yfinance as yf
from datetime import datetime, timedelta

def fetch_stock_data(symbol, period='1y'):
    """
    Fetch stock data for a given symbol from Yahoo Finance
    Returns data in format: date, stock, price
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

def fetch_multiple_stocks(symbols, period='1y'):
    """
    Fetch data for multiple stocks and combine into one DataFrame
    """
    all_data = []
    
    for symbol in symbols:
        print(f"Fetching data for {symbol}...")
        df = fetch_stock_data(symbol, period)
        if df is not None:
            all_data.append(df)
    
    if all_data:
        combined_df = pd.concat(all_data, ignore_index=True)
        return combined_df.sort_values(['date', 'stock'])
    
    return pd.DataFrame()

if __name__ == "__main__":
    # Example: Fetch data for a few popular S&P 500 stocks
    example_stocks = ['AAPL', 'MSFT', 'GOOGL', 'AMZN', 'TSLA']
    
    print("Fetching stock data for example stocks...")
    df = fetch_multiple_stocks(example_stocks, period='1mo')  # Last month for faster demo
    
    if not df.empty:
        print(f"\nFetched {len(df)} data points for {df['stock'].nunique()} stocks")
        print("\nFirst 10 rows of data:")
        print(df.head(10))
        
        print("\nData structure:")
        print(df.info())
        
        print("\nSample data by stock:")
        for stock in df['stock'].unique():
            stock_data = df[df['stock'] == stock]
            latest_price = stock_data['price'].iloc[-1]
            earliest_price = stock_data['price'].iloc[0]
            change = latest_price - earliest_price
            change_pct = (change / earliest_price) * 100
            print(f"{stock}: ${latest_price:.2f} (change: ${change:.2f}, {change_pct:.2f}%)")
    else:
        print("No data fetched!")