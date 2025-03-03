using Sequence.Position;
using Sequence.Position.Extract;

namespace PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow
{
    /// <summary>
    /// SlidingWindow位置クリエーター
    /// </summary>
    internal class SlidingWindowPositionCreator
    {
        private static readonly GenomePositionComparer _genomePositionComparer = new();

        private readonly SlidingWindowAnalysisSettings _settings;

        /// <summary>
        /// SlidingWindow位置クリエーターを作成する。
        /// </summary>
        /// <param name="settings">設定</param>
        public SlidingWindowPositionCreator(SlidingWindowAnalysisSettings settings)
        {
            _settings = settings;
        }

        /// <summary>
        /// SlidingWindow位置を作成する。
        /// </summary>
        /// <param name="variants">変異情報</param>
        /// <returns>SlidingWindow位置リスト</returns>
        public GenomePosition[] Create(IHasGenomePositionItem[] variants)
        {
            var windowPositionList = new List<GenomePosition>();

            var chrLookupQuery = variants
                .Select(x => x.GenomePosition)
                .ToLookup(x => x.ChrName);

            foreach (var gr in chrLookupQuery)
            {
                var lastPosition = gr.OrderBy(x => x, _genomePositionComparer).Last();
                var positions = Create(lastPosition);
                windowPositionList.AddRange(positions);
            }

            return [.. windowPositionList.OrderBy(x => x, _genomePositionComparer)];
        }

        /// <summary>
        /// １染色体分のWindow位置を作成する。
        /// </summary>
        /// <param name="lastVariantPosition">最終変異の位置情報</param>
        /// <returns>SlidingWindow位置リスト</returns>
        private GenomePosition[] Create(GenomePosition lastVariantPosition)
        {
            var positionList = new List<GenomePosition>();

            var chrName = lastVariantPosition.ChrName;
            var windowSize = _settings.WindowSize.BpValue;
            var stepSize = _settings.StepSize.BpValue;
            var currentStart = 1;
            var currentEnd = 1;

            while (currentEnd < lastVariantPosition.End)
            {
                currentEnd = currentStart + windowSize - 1;
                var windowPosition = new GenomePosition(chrName, currentStart, currentEnd);
                positionList.Add(windowPosition);

                currentStart += stepSize;
            }

            return [ ..positionList];
        }
    }
}
