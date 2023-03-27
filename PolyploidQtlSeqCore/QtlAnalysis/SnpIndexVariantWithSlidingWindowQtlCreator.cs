using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis
{
    /// <summary>
    /// SnpIndexVariantWithSlidingWindowQtlクリエーター
    /// </summary>
    internal class SnpIndexVariantWithSlidingWindowQtlCreator
    {
        private readonly GenomePositionExtractor<Window> _windowExtractor;
        private readonly ThreadNumber _threadNumber;

        /// <summary>
        /// SnpIndexVariantWithSlidingWindowQtlクリエーターを作成する。
        /// </summary>
        /// <param name="windows">全てのWindow</param>
        /// <param name="threadNumber">スレッド数</param>
        public SnpIndexVariantWithSlidingWindowQtlCreator(Window[] windows, ThreadNumber threadNumber)
        {
            if (windows.Length == 0) throw new ArgumentException(null, nameof(windows));

            _windowExtractor = new GenomePositionExtractor<Window>(windows);
            _threadNumber = threadNumber;
        }

        /// <summary>
        /// SnpIndexVariantWithSlidingWindowQtlを作成する。
        /// </summary>
        /// <param name="variants">QTL情報を持つ変異</param>
        /// <returns>SnpIndexVariantWithSlidingWindowQtl配列</returns>
        public SnpIndexVariantWithSlidingWindowQtl[] Create(SnpIndexVariantWithQtl[] variants)
        {
            var pOption = new ParallelOptions() { MaxDegreeOfParallelism = _threadNumber.Value };
            var results = new SnpIndexVariantWithSlidingWindowQtl[variants.Length];

            Parallel.For(0, variants.Length, pOption, i =>
            {
                var variant = variants[i];
                var overlapWindows = _windowExtractor.ExtractOverlap(variant.GenomePosition);
                var windows = new Windows(overlapWindows);
                var maxScoreWindowQtl = windows.ToMaxScoreSlidingWindowQtl();

                results[i] = new SnpIndexVariantWithSlidingWindowQtl(variant, maxScoreWindowQtl);
            });

            return results;
        }
    }
}
