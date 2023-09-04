using PolyploidQtlSeqCore.QtlAnalysis;

namespace PolyploidQtlSeq.Options.QtlAnalysis
{
    /// <summary>
    /// QTL解析シナリオ設定値DTO
    /// </summary>
    internal class QtlAnalysisScenarioSettingValueDTO : IQtlAnalysisScenarioSettingValue
    {
        /// <summary>
        /// QTL解析シナリオ設定値DTOインスタンスを作成する。
        /// </summary>
        public QtlAnalysisScenarioSettingValueDTO()
        {
            OutputDir = "";
            DisplayAnnotationImpacts = "";
            ThreadNumber = 0;
            Parent1MostAlleleRateThreshold = 0;
            Parent2SnpIndexRange = "";
            MinimumDepthThreshold = 0;
            MaxBulkSnpIndexThreshold = 0;
            Ploidy = 0;
            Parent2PlexNumber = 0;
            Bulk1Number = 0;
            Bulk2Number = 0;
            ReplicatesNumber = 0;
            WindowSize = 0;
            StepSize = 0;
            FigureHeight = 0;
            FigureWidth = 0;
            XAxisMajorStep = 0;
        }

        public string OutputDir { get; init; }

        public string DisplayAnnotationImpacts { get; init; }

        public int ThreadNumber { get; init; }

        public double Parent1MostAlleleRateThreshold { get; init; }

        public string Parent2SnpIndexRange { get; init; }

        public int MinimumDepthThreshold { get; init; }

        public double MaxBulkSnpIndexThreshold { get; init; }

        public int Ploidy { get; init; }

        public int Parent2PlexNumber { get; init; }

        public int Bulk1Number { get; init; }

        public int Bulk2Number { get; init; }

        public int ReplicatesNumber { get; init; }

        public int WindowSize { get; init; }

        public int StepSize { get; init; }

        public int FigureWidth { get; init; }

        public int FigureHeight { get; init; }

        public int XAxisMajorStep { get; init; }
    }
}
