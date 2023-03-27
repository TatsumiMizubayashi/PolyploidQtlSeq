using OxyPlot;
using OxyPlot.Series;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// 1染色体の変異コレクション
    /// </summary>
    internal class ChrVariants
    {
        private static readonly double _markerSize = 3;
        private static readonly MarkerType _markerType = MarkerType.Circle;

        private readonly SnpIndexVariantWithSlidingWindowQtl[] _variants;

        /// <summary>
        /// 1染色体の変異コレクションを作成する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <param name="variants">変異情報</param>
        public ChrVariants(string chrName, SnpIndexVariantWithSlidingWindowQtl[] variants)
        {
            if (string.IsNullOrWhiteSpace(chrName)) throw new ArgumentException(null, nameof(chrName));

            ChrName = chrName;
            _variants = variants;
        }

        /// <summary>
        /// 染色体名を取得する。
        /// </summary>
        public string ChrName { get; }

        /// <summary>
        /// Bulk1SnpIndexのScatterSeriesを作成する。
        /// </summary>
        /// <returns>Bulk1SnpIndexのScatterSeries</returns>
        public ScatterSeries CreateBulk1SnpIndexScatterSeries()
        {
            var series = CreateNewSeries(ColorPalette.Bulk1SnpIndexColor);

            var points = _variants
                .Where(x => x.Bulk1.SnpIndex.Value > 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Start.ToMbp(), x.Bulk1.SnpIndex.Value))
                .ToArray();

            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// Bulk1ZeroSnpIndexのScatterSeriesを作成する。
        /// </summary>
        /// <returns>Bulk1ZeroSnpIndexのScatterSeries</returns>
        public ScatterSeries CreateBulk1ZeroSnpIndexScatterSeries()
        {
            var series = CreateNewSeries(ColorPalette.Bulk1ZeroSnpIndexColor);

            var points = _variants
                .Where(x => x.Bulk1.SnpIndex.Value == 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Start.ToMbp(), x.Bulk1.SnpIndex.Value))
                .ToArray();

            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// Bulk2SnpIndexのScatterSeriesを作成する。
        /// </summary>
        /// <returns>Bulk2SnpIndexのScatterSeries</returns>
        public ScatterSeries CreateBulk2SnpIndexScatterSeries()
        {
            var series = CreateNewSeries(ColorPalette.Bulk2SnpIndexColor);

            var points = _variants
                .Where(x => x.Bulk2.SnpIndex.Value > 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Start.ToMbp(), x.Bulk2.SnpIndex.Value))
                .ToArray();

            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// Bulk2ZeroSnpIndexのScatterSeriesを作成する。
        /// </summary>
        /// <returns>Bulk2ZeroSnpIndexのScatterSeries</returns>
        public ScatterSeries CreateBulk2ZeroSnpIndexScatterSeries()
        {
            var series = CreateNewSeries(ColorPalette.Bulk2ZeroSnpIndexColor);

            var points = _variants
                .Where(x => x.Bulk2.SnpIndex.Value == 0)
                .Select(x => new ScatterPoint(x.GenomePosition.Start.ToMbp(), x.Bulk2.SnpIndex.Value))
                .ToArray();

            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        /// <summary>
        /// ΔSnpIndexのScatterSeriesを作成する。
        /// </summary>
        /// <returns>ΔSnpIndexのScatterSeries</returns>
        public ScatterSeries CreateDeltaSnpIndexScatterSeries()
        {
            var series = CreateNewSeries(ColorPalette.DeltaSnpIndexColor);

            var points = _variants
                .Select(x => new ScatterPoint(x.GenomePosition.Start.ToMbp(), x.DeltaSnpIndex.Value))
                .ToArray();

            if (points.Length > 0) series.Points.AddRange(points);

            return series;
        }

        private static ScatterSeries CreateNewSeries(OxyColor color)
        {
            return new ScatterSeries()
            {
                MarkerFill = color,
                MarkerSize = _markerSize,
                MarkerType = _markerType,
            };
        }
    }
}
