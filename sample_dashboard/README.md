# S&P 500 Stock Dashboard

A web-based dashboard built with Plotly Dash that displays stock price data for S&P 500 companies. The dashboard fetches real-time data from Yahoo Finance and allows users to select stocks and view price charts over the last year.

## Features

- **Stock Selection**: Dropdown menu with all S&P 500 companies
- **Interactive Charts**: Line chart showing stock price over the last year
- **Data Summary**: Table with key metrics including current price, 52-week high/low, and year-to-date performance
- **Real-time Data**: Fetches current stock data from Yahoo Finance
- **Responsive Design**: Clean, modern interface

## Data Structure

The dashboard processes data in the requested format where each row contains:
- `date`: Trading date
- `stock`: Stock symbol (S&P 500 companies)
- `price`: Closing price

## Installation

1. Navigate to the project directory:
```bash
cd sample_dashboard
```

2. Create and activate a virtual environment:
```bash
python3 -m venv venv
source venv/bin/activate  # On Windows: venv\Scripts\activate
```

3. Install the required dependencies:
```bash
pip install -r requirements.txt
```

4. Run the dashboard:
```bash
python dashboard.py
```

5. Open your web browser and navigate to:
```
http://localhost:8050
```

## Usage

### Running the Dashboard

1. **Select a Stock**: Use the dropdown menu to choose from any S&P 500 company
2. **View Chart**: The line chart will automatically update to show the selected stock's price over the last year
3. **Review Metrics**: Check the summary table below the chart for key performance indicators

### Example Data Script

You can also run the example script to see how the data is structured:

```bash
python example_usage.py
```

This script demonstrates how to fetch stock data in the required format (date, stock, price) for multiple S&P 500 companies.

## Technical Details

- **Framework**: Plotly Dash for the web application
- **Data Source**: Yahoo Finance API via `yfinance` library
- **Stock List**: S&P 500 companies scraped from Wikipedia (with fallback list)
- **Time Period**: Last 1 year of trading data
- **Update Frequency**: Real-time when stock selection changes

## Dependencies

- `dash`: Web application framework
- `plotly`: Interactive charting library
- `pandas`: Data manipulation
- `yfinance`: Yahoo Finance data fetching
- `requests`: HTTP requests for web scraping
- `beautifulsoup4`: HTML parsing for S&P 500 symbols

## Error Handling

The application includes robust error handling:
- Fallback list of popular S&P 500 stocks if web scraping fails
- Graceful handling of missing stock data
- Clear error messages for data fetch failures