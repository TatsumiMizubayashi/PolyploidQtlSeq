using OxyPlot;
using OxyPlot.Series;
using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// 1染色体のWindowコレクション
    /// </summary>
    internal class ChrWindows
    {
        private static readonly int _markerSize = 3;
        private static readonly MarkerType _lineMarkerType = MarkerType.None;
        private static readonly MarkerType _scatterMarkerType = MarkerType.Circle;

        private readonly Window[] _windows;

        /// <summary>
        /// 1染色体のWindowコレクションを作成する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <param name="windows">window</param>
        public ChrWindows(string chrName, Window[] windows)
        {
            if (string.IsNullOrWhiteSpace(chrName)) throw new ArgumentException(null, nameof(chrName));

            ChrName = chrName;
            _windows = windows;
        }

        /// <summary>
        /// 染色体名を取得する。
        /// </summary>
        public string ChrName { get; }

        /// <summary>
        /// Bulk1 SnpIndex平均値のLineSeriesを作成する。
        /// </summary>
        /// <returns>Bulk1 SnpIndex平均値のLineSeries</returns>
        public LineSeries CreateBulk1AverageSnpIndexLineSeries()
        {
            var series = CreateNewLineSeries(ColorPalette.AverageBulk1SnpIndexColor);
            var points = _windows
                .Select(x => new DataPoint(x.GenomePosition.Center.ToMbp(), x.AverageBulk1SnpIndex.Value))
                .ToArray();
            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// Bulk2 SnpIndex平均値のLineSeriesを作成する。
        /// </summary>
        /// <returns>Bulk2 SnpIndex平均値のLineSeries</returns>
        public LineSeries CreateBulk2AverageSnpIndexLineSeries()
        {
            var series = CreateNewLineSeries(ColorPalette.AverageBulk2SnpIndexColor);
            var points = _windows
                .Select(x => new DataPoint(x.GenomePosition.Center.ToMbp(), x.AverageBulk2SnpIndex.Value))
                .ToArray();
            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// ΔSnpIndex平均値のLineSeriesを作成する。
        /// </summary>
        /// <returns>ΔSnpIndex平均値のLineSeries</returns>
        public LineSeries CreateAverageDeltaSnpIndexLineSeries()
        {
            var series = CreateNewLineSeries(ColorPalette.AverageDeltaSnpIndexColor);
            var points = _windows
                .Select(x => new DataPoint(x.GenomePosition.Center.ToMbp(), x.AverageDeltaSnpIndex.Value))
                .ToArray();
            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// P95しきい値のLineSeriesを作成する。
        /// </summary>
        /// <returns>(+値のP95しきい値のLineSeries, -値のP95しきい値のLineSeries)</returns>
        public (LineSeries plusSeries, LineSeries minusSeries) CreateP95ThresholdLineSeries()
        {
            var plusSeries = CreateNewLineSeries(ColorPalette.P95Color);
            var plusPoints = _windows
                .Select(x => new DataPoint(x.GenomePosition.Center.ToMbp(), x.P95Qtl.ThresholdDeltaSnpIndex.Value))
                .ToArray();
            if (plusPoints.Length > 0) plusSeries.Points.AddRange(plusPoints);

            var minusSeries = CreateNewLineSeries(ColorPalette.P95Color);
            var minusPoints = plusPoints.Select(p => new DataPoint(p.X, p.Y * -1)).ToArray();
            if (minusPoints.Length > 0) minusSeries.Points.AddRange(minusPoints);

            return (plusSeries, minusSeries);
        }

        /// <summary>
        /// P99しきい値のLineSeriesを作成する。
        /// </summary>
        /// <returns>(+値のP99しきい値のLineSeries, -値のP99しきい値のLineSeries)</returns>
        public (LineSeries plusSeries, LineSeries minusSeries) CreateP99ThresholdLineSeries()
        {
            var plusSeries = CreateNewLineSeries(ColorPalette.P99Color);
            var plusPoints = _windows
                .Select(x => new DataPoint(x.GenomePosition.Center.ToMbp(), x.P99Qtl.ThresholdDeltaSnpIndex.Value))
                .ToArray();
            if (plusPoints.Length > 0) plusSeries.Points.AddRange(plusPoints);

            var minusSeries = CreateNewLineSeries(ColorPalette.P99Color);
            var minusPoints = plusPoints.Select(p => new DataPoint(p.X, p.Y * -1)).ToArray();
            if (minusPoints.Length > 0) minusSeries.Points.AddRange(minusPoints);

            return (plusSeries, minusSeries);
        }

        /// <summary>
        /// スコア(-log10(p))のScatterSeriesを作成する。
        /// </summary>
        /// <returns>スコア(-log10(p))のScatterSeries</returns>
        public ScatterSeries CreateScorexScatterSeries()
        {
            var series = CreateNewScatterSeries(ColorPalette.ScoreColor);

            var points = _windows
                .Select(x => new ScatterPoint(x.GenomePosition.Center.ToMbp(), x.AverageScore.Value))
                .ToArray();
            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// P95QTL数のScatterSeriesを作成する。
        /// </summary>
        /// <returns>P95QTL数のScatterSeries</returns>
        public ScatterSeries CreateP95QtlCountScatterSeries()
        {
            var series = CreateNewScatterSeries(ColorPalette.P95Color);

            var plusPoints = _windows
                .Where(x => x.P95Qtl.PlusDeltaSnpIndexQtlVariantCount > 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Center.ToMbp(), x.P95Qtl.PlusDeltaSnpIndexQtlVariantCount))
                .ToArray();
            if (plusPoints.Length > 0) series.Points.AddRange(plusPoints);

            var minusPoints = _windows
                .Where(x => x.P95Qtl.MinusDeltaSnpIndexQtlVariantCount > 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Center.ToMbp(), x.P95Qtl.MinusDeltaSnpIndexQtlVariantCount * -1))
                .ToArray();
            if (minusPoints.Length > 0) series.Points.AddRange(minusPoints);

            return series;
        }

        /// <summary>
        /// P99QTL数のScatterSeriesを作成する。
        /// </summary>
        /// <returns>P99QTL数のScatterSeries</returns>
        public ScatterSeries CreateP99QtlCountScatterSeries()
        {
            var series = CreateNewScatterSeries(ColorPalette.P99Color);

            var plusPoints = _windows
                .Where(x => x.P99Qtl.PlusDeltaSnpIndexQtlVariantCount > 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Center.ToMbp(), x.P99Qtl.PlusDeltaSnpIndexQtlVariantCount))
                .ToArray();
            if (plusPoints.Length > 0) series.Points.AddRange(plusPoints);

            var minusPoints = _windows
                .Where(x => x.P99Qtl.MinusDeltaSnpIndexQtlVariantCount > 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Center.ToMbp(), x.P99Qtl.MinusDeltaSnpIndexQtlVariantCount * -1))
                .ToArray();
            if (minusPoints.Length > 0) series.Points.AddRange(minusPoints);

            return series;
        }

        private static LineSeries CreateNewLineSeries(OxyColor color)
        {
            return new LineSeries()
            {
                MarkerFill = color,
                MarkerSize = _markerSize,
                MarkerType = _lineMarkerType,
                Color = color
            };
        }

        private static ScatterSeries CreateNewScatterSeries(OxyColor color)
        {
            return new ScatterSeries()
            {
                MarkerFill = color,
                MarkerSize = _markerSize,
                MarkerType = _scatterMarkerType
            };
        }
    }
}
