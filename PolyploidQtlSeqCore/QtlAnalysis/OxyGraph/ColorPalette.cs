using OxyPlot;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// カラーパレット
    /// </summary>
    internal static class ColorPalette
    {
        /// <summary>
        /// コンストラクタ
        /// </summary>
        static ColorPalette()
        {
            Bulk1SnpIndexColor = OxyColors.LightGreen;
            Bulk1ZeroSnpIndexColor = OxyColors.Green;
            AverageBulk1SnpIndexColor = OxyColors.Green;

            Bulk2SnpIndexColor = OxyColor.FromRgb(255, 184, 51);
            Bulk2ZeroSnpIndexColor = OxyColors.OrangeRed;
            AverageBulk2SnpIndexColor = OxyColors.OrangeRed;

            DeltaSnpIndexColor = OxyColor.FromRgb(187, 200, 230);
            AverageDeltaSnpIndexColor = OxyColor.FromRgb(25, 68, 142);

            ScoreColor = OxyColor.FromRgb(25, 68, 142);

            P95Color = OxyColor.FromRgb(248, 181, 0);
            P99Color = OxyColor.FromRgb(201, 23, 30);

            MajorGridlineColor = OxyColors.LightGray;
            GraphBackgroundColor = OxyColors.White;
        }

        /// <summary>
        /// Bulk1 SnpIndexの色
        /// </summary>
        public static OxyColor Bulk1SnpIndexColor { get; }

        /// <summary>
        /// Bulk1 SnpIndex=0の色
        /// </summary>
        public static OxyColor Bulk1ZeroSnpIndexColor { get; }

        /// <summary>
        /// Bulk1 SnpIndex平均値の色
        /// </summary>
        public static OxyColor AverageBulk1SnpIndexColor { get; }

        /// <summary>
        /// Bul2 SnpIndexの色
        /// </summary>
        public static OxyColor Bulk2SnpIndexColor { get; }

        /// <summary>
        /// Bulk2 SnpIndex=0の色
        /// </summary>
        public static OxyColor Bulk2ZeroSnpIndexColor { get; }

        /// <summary>
        /// Bulk2 SnpIndex平均値の色
        /// </summary>
        public static OxyColor AverageBulk2SnpIndexColor { get; }

        /// <summary>
        /// ΔSnpIndexの色
        /// </summary>
        public static OxyColor DeltaSnpIndexColor { get; }

        /// <summary>
        /// ΔSnpIndex平均値の色
        /// </summary>
        public static OxyColor AverageDeltaSnpIndexColor { get; }

        /// <summary>
        /// スコア(-log10(p))の色
        /// </summary>
        public static OxyColor ScoreColor { get; }

        /// <summary>
        /// P95の色
        /// </summary>
        public static OxyColor P95Color { get; }

        /// <summary>
        /// P99の色
        /// </summary>
        public static OxyColor P99Color { get; }

        /// <summary>
        /// Grid線の色
        /// </summary>
        public static OxyColor MajorGridlineColor { get; }

        /// <summary>
        /// グラフの背景色
        /// </summary>
        public static OxyColor GraphBackgroundColor { get; }
    }
}
