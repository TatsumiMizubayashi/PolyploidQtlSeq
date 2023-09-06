using PolyploidQtlSeqCore.Share;
using Sequence.Position;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// スライディングウインドウアナライザー
    /// </summary>
    internal class SlidingWindowAnalyzer
    {
        private readonly SlidingWindowAnalysisSettings _settings;
        private readonly ThreadNumber _threadNumber;

        /// <summary>
        /// スライディングウインドウアナライザーを作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        /// <param name="threadNumber">スレッド数</param>
        public SlidingWindowAnalyzer(SlidingWindowAnalysisSettings settings, ThreadNumber threadNumber)
        {
            _settings = settings;
            _threadNumber = threadNumber;
        }

        /// <summary>
        /// スライディングウインドウ解析を行う。
        /// </summary>
        /// <param name="variants">変異情報</param>
        /// <returns>Window配列</returns>
        public Window[] Analyze(SnpIndexVariantWithQtl[] variants)
        {
            var windowPositions = CreateWindowPositons(variants);
            var windows = new Window[windowPositions.Length];

            var pOption = new ParallelOptions() { MaxDegreeOfParallelism = _threadNumber.Value };
            var variantExtractor = new GenomePositionExtractor<SnpIndexVariantWithQtl>(variants);

            Parallel.For(0, windows.Length, pOption, i =>
            {
                var windowPosition = windowPositions[i];
                var extractVariants = variantExtractor.ExtractOverlap(windowPosition);
                var variantsInWindow = new VariantsInWindow(windowPosition, extractVariants);
                windows[i] = variantsInWindow.ToWindow();
            });

            return windows;
        }

        /// <summary>
        /// SlidingWindow位置を作成する。
        /// window位置は並べ替え済み。
        /// </summary>
        /// <param name="variants">変異</param>
        /// <returns>window位置</returns>
        private GenomePosition[] CreateWindowPositons(SnpIndexVariantWithQtl[] variants)
        {
            var creator = new SlidingWindowPositionCreator(_settings);
            return creator.Create(variants);
        }
    }
}
