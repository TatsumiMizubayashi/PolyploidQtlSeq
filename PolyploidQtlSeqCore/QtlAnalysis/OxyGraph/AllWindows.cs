using PolyploidQtlSeqCore.QtlAnalysis.SlidingWindow;
using System.Collections.Frozen;

namespace PolyploidQtlSeqCore.QtlAnalysis.OxyGraph
{
    /// <summary>
    /// AllWindow
    /// </summary>
    internal class AllWindows
    {
        private readonly Window[] _windows;
        private readonly FrozenDictionary<string, Window[]> _chrWindowsDictionary;

        /// <summary>
        /// AllWindowを作成する。
        /// </summary>
        /// <param name="windows">window</param>
        public AllWindows(Window[] windows)
        {
            _windows = windows;
            _chrWindowsDictionary = _windows.ToLookup(x => x.GenomePosition.ChrName)
                .ToDictionary(x => x.Key, x => x.ToArray())
                .ToFrozenDictionary();
        }

        /// <summary>
        /// X軸設定を作成する。
        /// </summary>
        /// <param name="majorStep">X軸MajorStep(MB)</param>
        /// <returns>X軸範囲</returns>
        public XAxisConfig CreateXAsisConfig(XAxisMajorStep majorStep)
        {
            return XAxisConfigCreator.Create(_windows, majorStep);
        }

        /// <summary>
        /// スコアY軸設定を作成する。
        /// </summary>
        /// <returns>スコアY軸設定</returns>
        public YAxisConfig CreateWindowScoreYAxisConfig()
        {
            return WindowScoreYAxisConfigCreator.Create(_windows);
        }

        /// <summary>
        /// QTL数Y軸設定を作成する。
        /// </summary>
        /// <returns>QTL数Y軸設定</returns>
        public YAxisConfig CreateQtlCountYAxisConfig()
        {
            return QtlCountYAxisConfigCreator.Create(_windows);
        }

        /// <summary>
        /// 変異数0Windowを除いたChrWindowsを取得する。
        /// </summary>
        /// <param name="chrName">染色多名</param>
        /// <returns>ChrWindows</returns>
        public ChrWindows GetZeroExcludeWindows(string chrName)
        {
            var windows = _chrWindowsDictionary[chrName]
                .Where(x => x.VariantCount.Count > 0)
                .ToArray();

            return new ChrWindows(chrName, windows);
        }

        /// <summary>
        /// 変異数0で分割したChrWindowsリストを取得する。
        /// </summary>
        /// <param name="chrName">染色体名</param>
        /// <returns>ChrWindows</returns>
        public IReadOnlyList<ChrWindows> GetZeroSplitWindowsList(string chrName)
        {
            var windows = _chrWindowsDictionary[chrName];
            var splitWindowsList = Split0Window(windows);

            return [.. splitWindowsList.Select(x => new ChrWindows(chrName, x))];
        }

        /// <summary>
        /// 個体数0Windowの位置で配列を分割する。
        /// </summary>
        /// <param name="windows">windows</param>
        /// <returns>window配列リスト</returns>
        private static List<Window[]> Split0Window(Window[] windows)
        {
            var resultWindowsList = new List<Window[]>();

            var residualWindows = new List<Window>(windows);
            while (residualWindows.Count > 0)
            {
                var targetWindows = residualWindows.TakeWhile(x => x.VariantCount.Count > 0).ToArray();
                if (targetWindows.Length > 0)
                {
                    resultWindowsList.Add(targetWindows);
                    residualWindows.RemoveRange(0, targetWindows.Length);
                }

                var zeroWindows = residualWindows.TakeWhile(x => x.VariantCount.Count == 0).ToArray();
                if (zeroWindows.Length != 0) residualWindows.RemoveRange(0, zeroWindows.Length);
            }

            return resultWindowsList;
        }
    }
}
